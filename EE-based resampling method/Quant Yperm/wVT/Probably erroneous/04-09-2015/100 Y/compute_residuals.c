#include <stdio.h>
#include <stdlib.h>

#include "svdcomp.h"
#include "nrutil.h"
#include "math.h"
#include "brent.h"

#include "compute_residuals.h"

#define MISSVAL_FUNC_XI -999999

extern void dposv_(char *uplo, int *n, int *nrhs, double *A, int *lda,
	double *B, int *ldb, int *info);

void compute_phi_svd(struct FAMILY family, struct PHI_SVD phi_svd){
	int i, j;
	int n_person = family.n_person;
	double **phi = family.phi;
	double **omega = dmatrix(1, n_person, 1, n_person);

	for (i = 1; i <= n_person; i++) {
		for (j = 1; j <= n_person; j++) {
			if (i == j) omega[i][j] = phi[i][j] + 1.0;
			else omega[i][j] = 2.0 * phi[i][j];
		}
	}

	svdcomp(omega, n_person, n_person, phi_svd.lambda, phi_svd.u);
	free_dmatrix(omega, 1, n_person, 1, n_person);
}

double compute_gee_res(struct FAMILY *family, struct PHI_SVD phi_svd, struct DATA_STRUCT data_struct){
	int i, j, k, fam;
	int n_total = data_struct.n_total;
	int n_person = family[1].n_person;
	int n_cov = data_struct.n_cov;
	int n_fam = data_struct.n_fam;

	double *trait = family[1].trait;
	double **covars = family[1].cov;
	double *mu = family[1].mu_hat;

	double *lambda = phi_svd.lambda;
	double **u = phi_svd.u;
	double **sig = phi_svd.sigma;

	double *optim_betahat = dvector(1, n_cov);
	double optim_sigma_t_sq = -1.0;
	double optim_xi = -1.0; // Will store the optimal heritability
	double *beta_init = dvector(1, n_cov); // Stores initial estimates for beta_hat

	// Getting yaffec and Cov matrix for all indivudals in the sample
	double *yaffec = dvector(1, n_total);
	double **Cov = dmatrix(1, n_total, 1, n_cov);

	i = 1;
	for (fam = 1; fam <= n_fam; fam++){
		trait = family[fam].trait;
		covars = family[fam].cov;
		for (j = 1; j <= n_person; j++){
			yaffec[i] = trait[j];
			for (k = 1; k <= n_cov; k++) Cov[i][k] = covars[j][k];
			i++;
		}
	}

	// yaffec is binary trait for all indiv and Cov is covariate matrix (n_total x n_cov)
	compute_logistic_coeff(yaffec, Cov, n_total, n_cov, beta_init);

	estimate_gee_parameters(n_total, u, lambda, beta_init,
		family, n_cov, n_fam, &optim_xi, &optim_sigma_t_sq, optim_betahat);

	/*for (i = 1; i <= n_cov; i++) printf("%f\n", optim_betahat[i]);
	printf("%f\n", optim_xi);
	exit(1);
*/
	// Compute estimated disease probabilities for each family & its members
	for (fam = 1; fam <= n_fam; fam++) {
		covars = family[fam].cov;
		mu = family[fam].mu_hat;
		for (j = 1; j <= n_person; j++) {
			mu[j] = 0.0;
			for (k = 1; k <= n_cov; k++) {
				mu[j] += optim_betahat[k] * covars[j][k];
			}
			mu[j] = exp(mu[j]) / (1 + exp(mu[j]));
		}
	}

	//compute Sigma and store it in phi_svd structure
	double **phi = family[1].phi;
	for (i = 1; i <= n_person; i++) {
		for (j = 1; j <= n_person; j++) {
			if (i == j) sig[i][j] = optim_xi*(phi[i][j] + 1.0) + (1.0 - optim_xi);
			else sig[i][j] = optim_xi*(2.0 * phi[i][j]);
		}
	}

	free_dvector(optim_betahat, 1, n_cov);
	free_dvector(beta_init, 1, n_cov);
	free_dvector(yaffec, 1, n_total);
	free_dmatrix(Cov, 1, n_total, 1, n_cov);

	return optim_xi;
}

// Computes estimated covariate effects assuming Yi's are independent (standard logistic model)
void compute_logistic_coeff(double *yaffec, double **Cov, int totpc, int n_cov, double *beta_init) {
	int i, j, k, l;

	int n = n_cov, nrhs = n_cov, lda = n_cov, ldb = n_cov, info;
	char *uplo = "L";
	double *a = dvector(0, n_cov*n_cov);
	double *b = dvector(0, n_cov*n_cov);

	double* pi = dvector(1, totpc);
	double* piold = dvector(1, totpc);
	double* w = dvector(1, totpc);
	double* z = dvector(1, totpc);
	for (i = 1; i <= totpc; i++){
		pi[i] = (yaffec[i] + 1.0 / 2.0) / 2.0;
		piold[i] = 0.0;
		w[i] = pi[i] * (1 - pi[i]);
		z[i] = (yaffec[i] - pi[i]) / w[i] + log(pi[i] / (1 - pi[i]));
	}

	double** XtWX = dmatrix(1, n_cov, 1, n_cov);  // Here X is D_\beta = \Gamma X
	double**XtWXinver = dmatrix(1, n_cov, 1, n_cov); // = D_\beta'\Gamma^-1 D_\beta = X'\Gamma X
	double* XtWZ = dvector(1, n_cov);

	double* betaLold = dvector(1, n_cov);
	double* betaLnew = dvector(1, n_cov);

	double diff = 1.0;
	double diffpi = 1.0;


	while (diffpi > 1e-6) {

		for (i = 1; i <= n_cov; i++)
			betaLold[i] = betaLnew[i];

		for (i = 1; i <= n_cov; i++) {
			for (j = 1; j <= n_cov; j++) {
				XtWX[i][j] = 0;
				for (k = 1; k <= totpc; k++) {
					XtWX[i][j] += Cov[k][i] * w[k] * Cov[k][j];
				}
			}
		}

		/*printf("\n%f \n",XtWX[1][1]);
		exit(1);*/

		for (j = 1; j <= n_cov; j++) {
			for (i = 1; i <= n_cov; i++) {
				a[(j - 1) * n_cov + (i - 1)] = XtWX[i][j]; // Indices go from 0 to 8
			}
		}

		for (j = 1; j <= n_cov; j++) {
			for (i = 1; i <= n_cov; i++) {
				b[(j - 1) * n_cov + (i - 1)] = ((i == j) ? 1 : 0);
			}
		}

		n = n_cov, nrhs = n_cov, lda = n_cov, ldb = n_cov;

		dposv_(uplo, &n, &nrhs, a, &lda, b, &ldb, &info);

		for (j = 1; j <= n_cov; j++) {
			for (i = 1; i <= n_cov; i++) {
				XtWXinver[i][j] = b[(j - 1) * n_cov + (i - 1)];
			}
		}
		//   printf("\n%f ",XtWXinver[1][1]);
		//    exit(1);
		for (i = 1; i <= n_cov; i++) {
			XtWZ[i] = 0;
			for (j = 1; j <= totpc; j++) {
				XtWZ[i] += Cov[j][i] * w[j] * z[j];
			}
		}
		//   printf("\n%f ",XtWZ[1]);
		// exit(1);
		for (i = 1; i <= n_cov; i++) {
			betaLnew[i] = 0;
			for (j = 1; j <= n_cov; j++) {
				betaLnew[i] += XtWXinver[i][j] * XtWZ[j];
			}
		}
		//printf("\nbetaLnew is:%f ", betaLnew[1]);
		//exit(1);

		diff = 0.0;
		for (i = 1; i <= n_cov; i++) diff += fabs(betaLnew[i] - betaLold[i]);
		diffpi = 0.0;
		for (i = 1; i <= totpc; i++) diffpi += fabs(pi[i] - piold[i]);
		//printf("diff is:%f %f %f %f.\n",diff,diffpi,betaLnew[2] - betaLold[2],betaLnew[2]);
		if (diffpi<1e-6) {
			for (i = 1; i <= n_cov; i++) betaLnew[i] = betaLold[i];
			if (diff>1e-5) printf("Fail to convergence; fitted probabilities numerically 0 or 1 occurred!\n");
		}
		else {
			for (i = 1; i <= totpc; i++) {
				piold[i] = pi[i];
				pi[i] = 0;
				for (j = 1; j <= n_cov; j++) {
					pi[i] += Cov[i][j] * betaLnew[j];
				}
				pi[i] = exp(pi[i]) / (1 + exp(pi[i]));
			}
			//for (i=1;i<=totpc;i++) printf("%f ",pi[i]);
			//exit(1);
			for (i = 1; i <= totpc; i++)
				w[i] = pi[i] * (1 - pi[i]);
			for (i = 1; i <= totpc; i++)
				z[i] = (yaffec[i] - pi[i]) / w[i] + log(pi[i] / (1 - pi[i])); // Like \Gamma^(-1)(Y-mu_hat)+ X\beta_hat
			//for (i=1;i<=totpc;i++) printf("%f ",z[i]);
		}
	}

	/*for (i = 1; i <= n_cov; i++) {
		logistic_coeff_var[i] = XtWXinver[i][i];
		}
		*/  // No need to compute variance for beta_hat??

	for (i = 1; i <= n_cov; i++)
		betaLold[i] = betaLnew[i];
	for (i = 1; i <= n_cov; i++){
		beta_init[i] = betaLnew[i];
	}

	free_dvector(a, 0, n_cov*n_cov);
	free_dvector(b, 0, n_cov*n_cov);
	free_dvector(pi, 1, totpc);
	free_dvector(piold, 1, totpc);
	free_dvector(w, 1, totpc);
	free_dvector(z, 1, totpc);
	free_dmatrix(XtWX, 1, n_cov, 1, n_cov);
	free_dmatrix(XtWXinver, 1, n_cov, 1, n_cov);
	free_dvector(XtWZ, 1, n_cov);
	free_dvector(betaLold, 1, n_cov);
	free_dvector(betaLnew, 1, n_cov);
}

void estimate_gee_parameters(int totpc, double **V, double *S, double *beta_init,
struct FAMILY *family, int n_cov, int n_fam, double *xi, double *sigma_t_sq, double *betahat) {

	int i;
	double f_root = -1.0;
	double xi_lo = 0, f_lo;
	double *current_beta = dvector(1, n_cov);
	double current_sigma_t_sq;
	param_func_xi parameters = { totpc, V, S, beta_init, family,
		current_beta, &current_sigma_t_sq, n_fam, n_cov };
	function_xi *func_xi = malloc(sizeof(function_xi));
	func_xi->parameters = parameters;
	func_xi->func_of_xi = &eqn3_lhs_minus_rhs;

	f_lo = eqn3_lhs_minus_rhs(0, parameters);

	if (f_lo <= 0.0) {
		//printf("Can't bracket the root. Estimate for xi is 0.0!");
		*xi = 0;
		*sigma_t_sq = current_sigma_t_sq;
		for (i = 1; i <= n_cov; i++) betahat[i] = current_beta[i];;
	}
	else {
		brent(func_xi, xi);
		//*xi=0.449952;
		f_root = eqn3_lhs_minus_rhs(*xi, parameters);
		*sigma_t_sq = current_sigma_t_sq;
		//printf("\ncurrent value is %f.", current_sigma_t_sq);
		for (i = 1; i <= n_cov; i++) betahat[i] = current_beta[i];
	}

	free_dvector(current_beta, 1, n_cov);
	free(func_xi);
}

double eqn3_lhs_minus_rhs(double xi, param_func_xi parameters) {
	return target_func_xi(xi, parameters.totpc, parameters.V, parameters.S,
		parameters.beta_init, parameters.family, parameters.current_beta,
		parameters.current_sigma_t_sq, parameters.n_fam, parameters.n_cov);
}

// totpc is n_total, xi is heritability
double target_func_xi(double xi, int totpc, double **V, double *S, double *beta_init,
struct FAMILY *family, double *current_beta, double *current_sigma_t_sq, int n_fam, int n_cov) {

	int i, j, k, l, fam;
	double *betahat;
	double sigma_t_sq;
	double rhs_eq3, lhs_eq3;
	double fvalue;
	int n_person = family[1].n_person;
	double **cov = family[1].cov;
	double **phi = family[1].phi;
	double *trait = family[1].trait;

	double **mu_f;
	mu_f = (double **)malloc((size_t)((n_fam + 1)*sizeof(double *)));
	mu_f[0] = NULL;

	for (fam = 1; fam <= n_fam; fam++) mu_f[fam] = dvector(1, n_person);

	//compute Sigma_inv; -> SAME FOR ALL FAMILIES!!
	double **Sigma_inv = dmatrix(1, n_person, 1, n_person);
	for (j = 1; j <= n_person; j++) {
		for (k = 1; k <= n_person; k++) {
			Sigma_inv[j][k] = 0.0;
			for (l = 1; l <= n_person; l++) {
				Sigma_inv[j][k] = Sigma_inv[j][k] + V[j][l] / (xi*S[l] + 1.0 - xi)*V[k][l];
			}
		}
	}
	//printf("\nNow We check Sigma_inv_f:\n"); // Should be indentity when xi=0
	//for (j = 1; j <= 6; j++) {
	//	for (k = 1; k <= 6; k++) {
	//		printf("%.2lf ", Sigma_inv_f[j][k]);
	//	}
	//	printf("\n");
	//}
	//exit(1);

	betahat = dvector(1, n_cov);
	int flag_betahat = -1;
	flag_betahat = compute_gee_betahat(Sigma_inv, beta_init, xi, betahat,
		family, mu_f, n_fam, n_cov);

	if (flag_betahat == 0) {
		fvalue = MISSVAL_FUNC_XI;
	}
	else {
		//compute mu_f
		for (fam = 1; fam <= n_fam; fam++) {
			cov = family[fam].cov;
			for (j = 1; j <= n_person; j++) {
				mu_f[fam][j] = 0.0;
				for (k = 1; k <= n_cov; k++) mu_f[fam][j] += betahat[k] * cov[j][k];
				mu_f[fam][j] = exp(mu_f[fam][j]) / (1 + exp(mu_f[fam][j]));
			}
		}
		//Assume sigma_t_sq (total variance)
		sigma_t_sq = 1.0;

		// Check to see whether these are the same in all families
		//Sigma_inverse * (Phi - I)
		double ***xi_V1 = (double ***)malloc((size_t)((n_fam + 1)*sizeof(double **)));
		xi_V1[0] = NULL;
		for (fam = 1; fam <= n_fam; fam++)
			xi_V1[fam] = dmatrix(1, n_person, 1, n_person);

		//Sigma_inverse * (Phi - I) * Sigma_inverse
		double ***xi_V2 = (double ***)malloc((size_t)((n_fam + 1)*sizeof(double **)));
		xi_V2[0] = NULL;
		for (fam = 1; fam <= n_fam; fam++)
			xi_V2[fam] = dmatrix(1, n_person, 1, n_person);


		//compute xi_V1 = Sigma_inverse * (Phi - I)
		for (fam = 1; fam <= n_fam; fam++) {
			phi = family[fam].phi;
			for (j = 1; j <= n_person; j++) {
				for (k = 1; k <= n_person; k++) {
					xi_V1[fam][j][k] = 0.0;
					for (l = 1; l <= n_person; l++) {
						if (l == k) xi_V1[fam][j][k] += Sigma_inv[j][l] * ((1 + phi[l][k]) - 1);
						else  xi_V1[fam][j][k] += Sigma_inv[j][l] * (2 * phi[l][k]);
					}
				}
			}
		}
		// Refers back to eqn 10 on sheng paper
		rhs_eq3 = 0.0;
		for (fam = 1; fam <= n_fam; fam++) {
			for (j = 1; j <= n_person; j++) {
				rhs_eq3 += xi_V1[fam][j][j];
			}
		}

		//compute xi_V2= Sigma_inverse * (Phi - I) * Sigma_inverse
		for (fam = 1; fam <= n_fam; fam++) {
			for (j = 1; j <= n_person; j++) {
				for (k = 1; k <= n_person; k++) {
					xi_V2[fam][j][k] = 0.0;
					for (l = 1; l <= n_person; l++) {
						xi_V2[fam][j][k] += xi_V1[fam][j][l] * Sigma_inv[l][k];
					}
				}
			}
		}
		//compute lhs_eq3
		lhs_eq3 = 0.0;
		for (fam = 1; fam <= n_fam; fam++) {
			trait = family[fam].trait;
			for (j = 1; j <= n_person; j++) {
				for (k = 1; k <= n_person; k++) {
					lhs_eq3 += (trait[j] - mu_f[fam][j]) / sqrt(mu_f[fam][j] * (1 - mu_f[fam][j])*mu_f[fam][k] * (1 - mu_f[fam][k]))*xi_V2[fam][j][k] * (trait[k] - mu_f[fam][k]);
				}
			}
		}
		//printf("%lf\n", lhs_eq3);

		//lhs_eq3 = lhs_eq3 / sigma_t_sq;
		fvalue = lhs_eq3 - rhs_eq3;
		for (i = 1; i <= n_cov; i++) current_beta[i] = betahat[i];
		*current_sigma_t_sq = sigma_t_sq;

		for (fam = 1; fam <= n_fam; fam++) {
			free_dmatrix(xi_V1[fam], 1, n_person, 1, n_person);
			free_dmatrix(xi_V2[fam], 1, n_person, 1, n_person);
		}
		free(xi_V1);
		free(xi_V2);
	}

	for (fam = 1; fam <= n_fam; fam++) free_dvector(mu_f[fam], 1, n_person);

	free_dmatrix(Sigma_inv, 1, n_person, 1, n_person);
	free(mu_f);
	free_dvector(betahat, 1, n_cov);

	return fvalue;
}

int compute_gee_betahat(double **Sigma_inv, double *beta_init, double xi, double *betahat, struct FAMILY *family, double **mu_f, int n_fam, int n_cov) {

	int i, j, k, l, fam, news;
	double diff_beta = 1.0;
	// printf("Current xi is: %f.\n",xi);
	int loop_while = 0;
	int flag_while = 1;

	double **cov = family[1].cov;
	double *trait = family[1].trait;
	int n_person = family[1].n_person;

	double *betaLold = dvector(1, n_cov);
	double *betaLnew = dvector(1, n_cov);
	// Let V= \Gamma^1/2 \Sigma^-1 \Gamma^1/2 = Var(Y)
	// X represents D_\beta = \Gamma X
	double **inv_XVinverX = dmatrix(1, n_cov, 1, n_cov); // (X'V^-1X)^-1 = sum_k (X_k'V^-1X_k) where k is family index
	double *XVinverY = dvector(1, n_cov); // X'V^-1Y= sum_k (X_k'V^-1Y_k) where k is family index

	double ***VinverX_f = (double ***)malloc((size_t)((n_fam + 1)*sizeof(double **)));
	VinverX_f[0] = NULL;
	for (fam = 1; fam <= n_fam; fam++)  VinverX_f[fam] = dmatrix(1, n_person, 1, n_cov);// V^-1X (n x k) matrix

	double *a;
	a = (double *)malloc((size_t)((n_cov * n_cov + 1) * sizeof(double)));
	double *b;
	b = (double *)malloc((size_t)((n_cov * n_cov + 1) * sizeof(double)));

	for (i = 1; i <= n_cov; i++) {
		betaLnew[i] = beta_init[i];
		betaLold[i] = beta_init[i];
	}

	while (diff_beta > 10e-9) {

		loop_while++;
		if (loop_while >= 1000) {
			FILE *optim_error;
			optim_error = fopen("optim_error.txt", "a+");
			fprintf(optim_error, "The while loop in function compute_gee_betahat"
				" is deadlocked, where xi is %f.\n", xi);
			fclose(optim_error);
			flag_while = 0;
			break;
		}
		for (i = 1; i <= n_cov; i++) betaLold[i] = betaLnew[i];

		//compute mu_f
		for (fam = 1; fam <= n_fam; fam++) {
			cov = family[fam].cov;
			for (j = 1; j <= n_person; j++) {
				mu_f[fam][j] = 0.0;
				for (k = 1; k <= n_cov; k++) mu_f[fam][j] += betaLold[k] * cov[j][k];
				mu_f[fam][j] = exp(mu_f[fam][j]) / (1 + exp(mu_f[fam][j]));
			}
		}
		/*printf("\ncheck mu_f:\n");
		i = 1;
		for (j = 1; j <= totpc_f[i]; j++) {
		printf("%f ", mu_f[i][j]);
		printf("\n");
		}*/

		// Compute V^-1X for each family
		for (fam = 1; fam <= n_fam; fam++) {
			cov = family[fam].cov;
			for (j = 1; j <= n_person; j++) {
				for (k = 1; k <= n_cov; k++) {
					VinverX_f[fam][j][k] = 0.0;
					for (l = 1; l <= n_person; l++){
						VinverX_f[fam][j][k] += Sigma_inv[j][l] * cov[l][k] * sqrt(mu_f[fam][l] * (1 - mu_f[fam][l])*mu_f[fam][j] * (1 - mu_f[fam][j]));
					}
				}
			}
		}
		/*for (i = 1; i <= 10; i++){
			printf("%f\t\t",mu_f[1][i]);
			printf("%d %f %f %f\n", i, VinverX_f[1][i][1], VinverX_f[1][i][2], VinverX_f[1][i][3]);
			}
			exit(1);
			*/
		for (i = 1; i <= n_cov; i++) {
			for (j = 1; j <= n_cov; j++) {
				inv_XVinverX[i][j] = 0.0;
				for (fam = 1; fam <= n_fam; fam++) {
					cov = family[fam].cov;
					for (l = 1; l <= n_person; l++) {
						inv_XVinverX[i][j] += cov[l][i] * VinverX_f[fam][l][j];
					}
				}
			}
		}
		//for (i = 1; i <= n_cov; i++){
		//	printf("%d %f %f %f\n", i, inv_XVinverX[i][1], inv_XVinverX[i][2], inv_XVinverX[i][3]);
		//}
		//exit(1);


		for (j = 1; j <= n_cov; j++) {
			for (i = 1; i <= n_cov; i++) {
				a[(j - 1) * n_cov + (i - 1)] = inv_XVinverX[i][j];
			}
		}

		for (j = 1; j <= n_cov; j++) {
			for (i = 1; i <= n_cov; i++) {
				b[(j - 1) * n_cov + (i - 1)] = (i == j) ? 1 : 0;
			}
		}

		int n = n_cov, nrhs = n_cov, lda = n_cov, ldb = n_cov, info;
		char *uplo = "L";
		dposv_(uplo, &n, &nrhs, a, &lda, b, &ldb, &info);

		for (j = 1; j <= n_cov; j++) {
			for (i = 1; i <= n_cov; i++) {
				inv_XVinverX[i][j] = b[(j - 1) * n_cov + (i - 1)];
			}
		}

		//compute XtVinverY
		for (fam = 1; fam <= n_fam; fam++) {
			cov = family[fam].cov;
			for (j = 1; j <= n_person; j++) {
				for (k = 1; k <= n_cov; k++) {
					VinverX_f[fam][j][k] = 0.0;
					for (l = 1; l <= n_person; l++){
						VinverX_f[fam][j][k] += Sigma_inv[j][l] * cov[l][k] / sqrt(mu_f[fam][j] * (1 - mu_f[fam][j]))*sqrt(mu_f[fam][l] * (1 - mu_f[fam][l]));
					}
				}
			}
		}

		//// This should print X for the first iteration (i.e. when xi=0)
		//for (i = 1; i <= n_person; i++) printf("%f %f %f\n", VinverX_f[1][i][1], VinverX_f[1][i][2], VinverX_f[1][i][3]);
		//exit(1);

		for (i = 1; i <= n_cov; i++){
			XVinverY[i] = 0.0;
			for (fam = 1; fam <= n_fam; fam++){
				trait = family[fam].trait;
				for (k = 1; k <= n_person; k++) {
					XVinverY[i] += VinverX_f[fam][k][i] * (trait[k] - mu_f[fam][k]);
				}
			}
		}
		/*for (i = 1; i <= n_cov; i++){ printf("%lf ", XVinverY[i]*10000000); }
		exit(1);*/

		//solve quasi-likelihood equations for beta using Newton's method with Fisher scoring
		for (i = 1; i <= n_cov; i++) {
			for (j = 1; j <= n_cov; j++) {
				betaLnew[i] += inv_XVinverX[i][j] * XVinverY[j];
			}
		}
		//printf("\ncheck betaLnew:");
		//for (i = 1; i <= n_cov; i++){ printf("%f ", betaLnew[i]); }
		//exit(1);

		//compute beta_diff
		diff_beta = 0.0;
		for (i = 1; i <= n_cov; i++) diff_beta += fabs(betaLnew[i] - betaLold[i]);
		//printf("\ndiff_beta is %f.",diff_beta);
		//  printf("Current loop is %d.\n",loop_while);
	}

	for (i = 1; i <= n_cov; i++) {
		if ((int)isnan(betaLnew[i]) == 1) flag_while = 0;
	}

	if (flag_while == 1) {
		for (i = 1; i <= n_cov; i++) betahat[i] = betaLnew[i];
	}
	else {
		for (i = 1; i <= n_cov; i++) betahat[i] = 0.0;
	}

	for (fam = 1; fam <= n_fam; fam++) free_dmatrix(VinverX_f[fam], 1, n_person, 1, n_cov);

	free_dvector(betaLold, 1, n_cov);
	free_dvector(betaLnew, 1, n_cov);
	free_dmatrix(inv_XVinverX, 1, n_cov, 1, n_cov);
	free(VinverX_f);
	free_dvector(XVinverY, 1, n_cov);
	free(a);
	free(b);

	return flag_while;
}




