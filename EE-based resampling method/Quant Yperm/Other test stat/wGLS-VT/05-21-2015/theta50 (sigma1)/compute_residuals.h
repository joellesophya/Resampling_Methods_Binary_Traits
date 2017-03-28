#ifndef COMPUTE_RESIDUALS_H
#define COMPUTE_RESIDUALS_H

#include "simulation.h"

void compute_phi_svd(struct FAMILY family, struct PHI_SVD phi_svd);
double compute_gee_res(struct FAMILY *family, struct PHI_SVD phi_svd, struct DATA_STRUCT data_struct);
void compute_logistic_coeff(double *yaffec, double **Cov, int totpc, int n_cov, double *beta_init);

void estimate_gee_parameters(int totpc, double **V, double *S, double *beta_init, struct FAMILY *family, int n_cov, int n_fam, double *xi, double *sigma_t_sq, double *betahat);

typedef struct {
	int totpc;
	double **V;
	double *S;
	double *beta_init;
	struct FAMILY *family;
	double *current_beta;
	double *current_sigma_t_sq;
	int n_fam;
	int n_cov;
} param_func_xi;
typedef struct {
	param_func_xi parameters;
	double(*func_of_xi)(double, param_func_xi);
} function_xi;
double eqn3_lhs_minus_rhs(double xi, param_func_xi parameters);
double target_func_xi(double xi, int totpc, double **V, double *S, double *beta_init, struct FAMILY *family, double *current_beta, double *current_sigma_t_sq, int n_fam, int n_cov);

int compute_gee_betahat(double **Sigma_inv, double *beta_init, double xi,
	double *betahat, struct FAMILY *family, double **mu_f, int n_fam, int n_cov);



#endif 
