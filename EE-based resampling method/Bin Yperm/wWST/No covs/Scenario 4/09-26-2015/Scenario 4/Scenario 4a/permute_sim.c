#include "permute_sim.h"

extern void dposv_(char *uplo, int *n, int *nrhs, double *A, int *lda,
	double *B, int *ldb, int *info);

void permute_pre(struct FAMILY *family, struct DATA_STRUCT data_struct, struct PHI_SVD phi_svd, struct TMP_PERM tmp_perm){
	int i, j, k, l, fam, tot;
	int n_fam = data_struct.n_fam;
	int n_total = data_struct.n_total;
	int n_cov = 1;
	int n_person = family[1].n_person;
	double tol = data_struct.tol;
	double tmp_hat = 0, tmp = 0;

	//1 Compute aug=[\Gamma^1/2X_f  \Gamma^{-1/2)(Y-mu_hat)] for each family
	// Get cholesky of omega in order to compute W=C{-t}aug
	// Also compute WtW
	double **sigma = phi_svd.sigma;
	double **chol = tmp_perm.cholesky;
	for (i = 1; i <= n_person; i++) {
		for (j = 1; j <= n_person; j++) tmp_perm.cholesky[i][j] = 0;
	}

	double **aug = dmatrix(1, n_person, 1, n_cov + 1);
	double **cholaug = dmatrix(1, n_person, 1, n_cov + 1);

	double *mu = family[1].mu_hat;
	double **cov = family[1].cov;
	double *trait = family[1].trait;

	double **WtW = dmatrix(1, n_cov, 1, n_cov);
	double **W = dmatrix(1, n_total, 1, n_cov);
	double *c_res = dvector(1, n_total);
	double **WtW_inv = dmatrix(1, n_cov, 1, n_cov);
	double *a = dvector(0, n_cov*n_cov);
	double *b = dvector(0, n_cov*n_cov);
	double **omega = dmatrix(1, n_total, 1, n_total);
	double **vs = dmatrix(1, n_total, 1, n_total);
	double *evals = dvector(1, n_total);



	// Initialize vectors
	for (i = 1; i <= n_total; i++){
		c_res[i] = 0;
		evals[i] = 0;
		for (j = 1; j <= n_cov; j++){
			W[i][j] = 0;
		}
		if (i == 1){
			for (j = 1; j <= n_cov*n_cov; j++){
				a[j] = 0;
				b[j] = 0;
			}
		}
	}

	for (i = 1; i <= n_cov; i++){
		for (j = 1; j <= n_cov; j++){
			WtW[i][j] = 0.0;
			WtW_inv[i][j] = 0.0;
			tot = 1;
			for (fam = 1; fam <= n_fam; fam++){
				mu = family[fam].mu_hat;
				cov = family[fam].cov;
				trait = family[fam].trait;

				// Compute aug for the family 
				for (k = 1; k <= n_person; k++) {
					for (l = 1; l <= n_cov; l++) {
						aug[k][l] = cov[k][l] * sqrt(mu[k] * (1 - mu[k]));
					}
					aug[k][n_cov + 1] = (trait[k] - mu[k]) / sqrt(mu[k] * (1 - mu[k]));
					if ((i == 1)&(j == 1)){
						for (l = 1; l <= n_cov; l++) chol[k][l] == 0;
					}
				}

				k = cholesky(sigma, n_person, aug, n_cov + 1, chol, cholaug, 1);
				if (k == 0){
					printf("Non-psd matrix!");
					exit(1);
				}

				for (k = 1; k <= n_person; k++){
					WtW[i][j] += cholaug[k][i] * cholaug[k][j];
					//if (j == 1) tmp += cholaug[k][i] * cholaug[k][n_cov + 1]; // Check that null estimates were computed correctly using score eqn
					if ((i == 1) && (j == 1)){
						c_res[tot] = cholaug[k][n_cov + 1];
						for (l = 1; l <= n_cov; l++)
							W[tot][l] = cholaug[k][l];
						tot++;
					}
				}
			}
		}
	}
	//printf("%f", tmp);
	//exit(1);

	//2 Compute (WtW)_inv 
	int n = n_cov, nrhs = n_cov, lda = n_cov, ldb = n_cov, info;
	char *uplo = "L";

	for (j = 1; j <= n_cov; j++) {
		for (i = 1; i <= n_cov; i++) {
			a[(j - 1) * n_cov + (i - 1)] = WtW[i][j]; // Indices go from 0 to n_cov-1
			b[(j - 1) * n_cov + (i - 1)] = ((i == j) ? 1 : 0);
		}
	}

	n = n_cov, nrhs = n_cov, lda = n_cov, ldb = n_cov;

	dposv_(uplo, &n, &nrhs, a, &lda, b, &ldb, &info);

	for (j = 1; j <= n_cov; j++) {
		for (i = 1; i <= n_cov; i++) {
			WtW_inv[i][j] = b[(j - 1) * n_cov + (i - 1)];
		}
	}

	//3 Get \Omega=I-W(WtW_inv)W'
	for (i = 1; i <= n_total; i++){
		for (j = 1; j <= n_total; j++){
			omega[i][j] = 0.0;
			vs[i][j] = 0;
			for (k = 1; k <= n_cov; k++){
				for (l = 1; l <= n_cov; l++){
					omega[i][j] += W[i][k] * WtW_inv[k][l] * W[j][l];
				}
			}
			if (i == j) omega[i][j] = 1.0 - omega[i][j];
			else omega[i][j] = (-1.0)*omega[i][j];
		}
	}

	//4 Compute V where \Omega= VDV' and V=[V1 V0] -> Take svd of \Omega
	svdcomp(omega, n_total, n_total, evals, vs);
	int nonzero = 0;
	for (i = 1; i <= n_total; i++){
		if (evals[i] > tol)
			nonzero++;
	}

	if (nonzero != (n_total - n_cov)){
		printf("Error in svd decomposition of Omega");
		exit(1);
	}

	//5 Obtain V1
	tot = 1;
	double **v1 = tmp_perm.v1;
	for (i = 1; i <= n_total; i++){
		if (evals[i] > tol){
			for (j = 1; j <= n_total; j++)	v1[j][tot] = vs[j][i];
			tot++;
		}
	}

	//6 Compute delta = V1'*c_res
	double *delta = tmp_perm.delta;
	for (i = 1; i <= nonzero; i++){
		delta[i] = 0.0;
		for (j = 1; j <= n_total; j++)
			delta[i] += v1[j][i] * c_res[j];
	}

	free_dmatrix(aug, 1, n_person, 1, n_cov + 1);
	free_dmatrix(cholaug, 1, n_person, 1, n_cov + 1);
	free_dvector(a, 0, n_cov*n_cov);
	free_dvector(b, 0, n_cov*n_cov);
	free_dmatrix(WtW, 1, n_cov, 1, n_cov);
	free_dmatrix(W, 1, n_total, 1, n_cov);
	free_dmatrix(WtW_inv, 1, n_cov, 1, n_cov);
	free_dmatrix(omega, 1, n_total, 1, n_total);
	free_dmatrix(vs, 1, n_total, 1, n_total);
	free_dvector(evals, 1, n_total);
	free_dvector(c_res, 1, n_total);
}

void permute(struct FAMILY *family, struct DATA_STRUCT data_struct, struct TMP_PERM tmp_perm) {
	int i, j, fam, ind_i, tot;

	int n_person = family[1].n_person; //Assume same number of indiv for all families
	int n_total = data_struct.n_total;
	int n_fam = data_struct.n_fam;
	int n_cov = 1;

	int nonzero = n_total - n_cov;

	double *delta = tmp_perm.delta;
	double **v1 = tmp_perm.v1;
	double **chol = tmp_perm.cholesky;
	double *y_perm = family[1].y_perm;
	double *mu = family[1].mu_hat;


	//7 Permute entries of delta
	shuffle(delta, nonzero);

	//8 Compute V1*delta

	for (i = 1; i <= n_total; i++){
		tmp_perm.v1Delta[i] = 0.0;
		for (j = 1; j <= nonzero; j++){
			tmp_perm.v1Delta[i] += v1[i][j] * delta[j];
		}
	}

	//9 Separate the entries of (V1*delta) into n_fam vectors of size n_person
	tot = 1;
	for (fam = 1; fam <= n_fam; fam++){
		for (i = 1; i <= n_person; i++){
			tmp_perm.v1Delta_f[fam][i] = tmp_perm.v1Delta[tot];
			tot++;
		}
	}

	//10 Multiply each (V1*delta)f by \Gamma^(1/2)Ct for each family
	double *Cp = dvector(1, n_total), thresh = 0, sum_y = 0;
	int g = 1;
	for (fam = 1; fam <= n_fam; fam++){
		y_perm = family[fam].y_perm;
		mu = family[fam].mu_hat;
		for (i = 1; i <= n_person; i++, g++){
			y_perm[i] = 0.0;
			for (j = 1; j <= n_person; j++) y_perm[i] += sqrt(mu[i] * (1 - mu[i]))*chol[j][i] * tmp_perm.v1Delta_f[fam][j];
			y_perm[i] += mu[i];
			
			Cp[g] = y_perm[i];
			sum_y += family[fam].trait[i];
		}
	}
	
	// double mu2_y = 0, tmp1 = 0, tmp_hat = 0.0;
	// for (fam = 1; fam <= n_fam; fam++){
		// for (i = 1; i <= n_person; i++){
			// tmp1 += family[fam].trait[i];
			// tmp_hat += family[fam].y_perm[i];
			// mu2_y += family[fam].y_perm[i] * family[fam].y_perm[i];
		// }
	// }

	// printf("%f %f\t", tmp_hat, tmp1);
	// printf("%f %f\n", mu2_y - tmp_hat*tmp_hat / n_total, tmp1 - (tmp1 * tmp1 / n_total));
	//exit(1);
	
	// Make y_perm binary
	qsort(&Cp[1], n_total, sizeof(double), cmpfunc);
	thresh = Cp[1 + (int) (n_total - sum_y)]; // Set top Y_bar replicates to 1, the rest to 0
	
	for (fam = 1; fam <= n_fam; fam++){
		for (i = 1; i <= n_person; i++){
			if (family[fam].y_perm[i]>= thresh) family[fam].y_perm[i] = 1.0;
			else family[fam].y_perm[i] = 0.0;
		}
	}
	
	free_dvector(Cp, 1, n_total);
}

int cmpfunc (const void * a, const void * b)
{
	if (*(double*)a > *(double*)b) return 1;
	else if (*(double*)a < *(double*)b) return -1;
	else return 0; 
}