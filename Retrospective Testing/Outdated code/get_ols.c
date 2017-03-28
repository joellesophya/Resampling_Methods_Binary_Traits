#include <stdio.h>
#include <stdlib.h>

#include "nrutil.h"
#include "math.h"
#include "cholesky.h"
#include "simulation.h"

#include "get_ols.h"

void get_ols_y(struct FAMILY *family, struct DATA_STRUCT data_struct){
	int i, j, k, j_cov;
	int n_person = family[1].n_person;
	int n_fam = data_struct.n_fam;
	int n_cov = data_struct.n_cov;
	int n_total = data_struct.n_total;
	
	double **xtxI = data_struct.xtxInv;
	double *sigma_y = data_struct.sig_hat;
	double **ww = dmatrix(1, n_cov, 1, n_cov);
	double *wd = dvector(1, n_cov);
	double *bhat = dvector(1, n_cov);
	
	for (i = 1; i <= n_cov; i++) {
		for (j = 1, wd[i] = 0; j <= n_cov; j++) ww[i][j] = 0;		
	}
	
	// Compute X^tX and X^tY
	for (i = 1; i <= n_fam; i++) {
		for (j = 1; j <= n_person; j++) {
			for (k = 1; k <= n_cov; k++) {
				wd[k] += family[i].cov[j][k] * family[i].trait[1][j];
				for (j_cov = 1; j_cov <= n_cov; j_cov++) ww[k][j_cov] += family[i].cov[j][k] * family[i].cov[j][j_cov];
			}
		}
	}
	
	double **chol = dmatrix(1, n_cov, 1, n_cov);
	double **aug = dmatrix(1, n_cov, 1, n_cov + 1);
	
	for (i = 1; i<= n_cov; i++) {
		for (j = 1; j <= n_cov; j++) {
			chol[i][j] = 0.0;
			aug[i][j] = 0.0;
		}
		aug[i][i] = 1.0;
		aug[i][n_cov + 1] = wd[i];
	}
	// Get (X^tX)^(-1)
	int posdef = 1;
	posdef = cholesky(ww, n_cov, aug, n_cov + 1, chol, aug, 1);
	if (posdef == 0) {
		printf("\nERROR: The matrix W^T W is not positive semi-definite\n"
			"(W is the matrix of covariates)\n");
		exit(1);
	}
	
	for (i = 1; i <= n_cov; i++) {		
		for (j = 1, bhat[i] = 0; j <= n_cov; j++) {
			bhat[i] += aug[j][i] * aug[j][n_cov + 1];
			for (k = 1, xtxI[i][j] = 0; k <= n_cov; k++) xtxI[i][j] += aug[k][i] * aug[k][j];
		}
	}
	
	// Store the OLS residuals and get MSE
	double tmp_e = 0;
	for (i = 1; i <= n_fam; i++) {		
		for (j = 1; j <= n_person; j++) {
			for (k = 1; k <= n_cov; k++) family[i].trait[1][j] -= family[i].cov[j][k] * bhat[k];
			tmp_e += family[i].trait[1][j] * family[i].trait[1][j];
		}
	}
	sigma_y[1] = tmp_e / (n_total - n_cov) ;
	
	free_dmatrix(ww, 1, n_cov, 1, n_cov);
	free_dvector(wd, 1, n_cov);
	free_dvector(bhat, 1, n_cov);
	free_dmatrix(chol, 1, n_cov, 1, n_cov);
	free_dmatrix(aug, 1, n_cov, 1, n_cov + 1);
	
}

void get_ols_y_rep(struct FAMILY *family, struct DATA_STRUCT data_struct, int i_perm){
	int i, j, k;
	int n_person = family[1].n_person;
	int n_fam = data_struct.n_fam;
	int n_cov = data_struct.n_cov;
	int n_total = data_struct.n_total;
	
	double **xtxI = data_struct.xtxInv;
	double *sigma_y = data_struct.sig_hat;
	double *wd = dvector(1, n_cov);
	double *bhat = dvector(1, n_cov);
	
	for (i = 1; i <= n_cov; i++) wd[i] = 0; 
	
	// Compute X^tY_pi
	for (i = 1; i <= n_fam; i++) {
		for (j = 1; j <= n_person; j++) {
			for (k = 1; k <= n_cov; k++) wd[k] += family[i].cov[j][k] * family[i].trait[i_perm + 1][j];
		}
	}
	
	// Get OLS b_hat
	for (i = 1; i <= n_cov; i++) {		
		for (j = 1, bhat[i] = 0; j <= n_cov; j++) bhat[i] += xtxI[i][j] * wd[j];			
	}
	
	// Store the OLS residuals and get MSE
	double tmp_e = 0;
	for (i = 1; i <= n_fam; i++) {		
		for (j = 1; j <= n_person; j++) {
			for (k = 1; k <= n_cov; k++) family[i].trait[i_perm + 1][j] -= family[i].cov[j][k] * bhat[k];
			tmp_e += family[i].trait[i_perm + 1][j] * family[i].trait[i_perm + 1][j];
		}
	}
	sigma_y[i_perm + 1] = tmp_e / (n_total - n_cov) ;
	
	free_dvector(wd, 1, n_cov);
	free_dvector(bhat, 1, n_cov);
}
