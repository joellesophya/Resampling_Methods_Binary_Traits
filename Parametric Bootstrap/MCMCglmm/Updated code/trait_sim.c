#include "trait_sim.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "nrutil.h"
#include "gadev.h"
#include "cholesky.h"
#include "rand1.h"

void run_trait_sim(struct FAMILY family, double *mles, int n_cov, int i_perm, long *seed) {
	int i, n_person = family.n_person;
	double *pi_mu = dvector(1, n_person);
	for (i = 1; i <= n_person; i++) pi_mu[i] = 0;

	// Sample random effects
	if (mles[n_cov + 1] > 0) mvn_ran_obs(pi_mu, n_person, family.phi, mles[n_cov + 1], seed);
	// For the fixed effects
	agesexeffect(pi_mu, family.cov, n_person, mles, n_cov); // Adds pi_mu= Xi'\beta 

	// Use logistic model
	compute_trait_logistic(family.trait[i_perm + 1], pi_mu, n_person, seed);

	free_dvector(pi_mu, 1, n_person);
}

//Subtracted one from sex variable to make it into an indicator var (0= male and 1=female)
// Added in age effect here as well
void agesexeffect(double *mu, double **cov, int n_person, double *mles, int n_cov) {
	int i, j;

	for (i = 1; i <= n_person; i++){
		for (j = 1; j <= n_cov; j++) mu[i] += cov[i][j] * mles [j];
	}
}

/* Program to simulate a multivariate normal distribution with dependence
   We want to simulate under a model multivariate normal model:

   MVN(fixed effects, Omega) = fixed effects + MVN(O, Omega)
   where Omega = Phi_RR*sigma_a^2

   Diagonal elements of Phi_RR are 1+h_i (h_i is the inbreeding coefficient)
   Off diagnoal elements are 2*phi_ij (phi_ij is the kinship coefficient)

   Since Omega is p.s.d. there exists a upper triangular matrix C s.t
   Omega = C^tC (Cholsky decomposition)

   Then we can simulate X as MVN(0,I) and get our desired Y vector
   as Y = C^tX ~ MVN(0,C^tIC) = MVN(0,Omega) and then add the mu:
   MVN(genetic effect, Omega) = mu + MVN(0, Omega) = MVN(mu, Omega)
   */
void mvn_ran_obs(double *mu, int n_person, double **phi, double add_var, long *seed) {
	int i, j;
	double *trait_tmp, **omega, **chol, tmp;

	trait_tmp = dvector(1, n_person);
	omega = dmatrix(1, n_person, 1, n_person);
	chol = dmatrix(1, n_person, 1, n_person);
	for (i = 1; i <= n_person; i++){
		for (j = 1; j <= n_person; j++) chol[i][j] = 0;
	}

	/* draw from N(0, I) */
	for (i = 1; i <= n_person; i++) trait_tmp[i] = gadev(seed);

	/* create Omega = sigma_a^2*Phi I*/
	omega_fun(omega, n_person, phi, add_var);

	/* find the cholesky for Omega = C^tC*/
	cholesky(omega, n_person, NULL, 0, chol, NULL, 1);

	/* Change the dependance structure by using the fact that if X~MVN(0,I),
	   then U = C^tX ~ MVN(0,C^t*I*C) = MVN(0, Omega)
	   For upper triangular matrix C, U = C^tX may be calculated as:
	   U_i = sum_{j=1 to j=i} C[j][i] * X_j	 */
	for (i = 1; i <= n_person; i++) {
		tmp = 0.0;
		for (j = 1; j <= i; j++) tmp = tmp + chol[j][i] * trait_tmp[j];
		mu[i] += tmp; /// Adds u_i to X_i'\beta
	}

	free_dvector(trait_tmp, 1, n_person);
	free_dmatrix(omega, 1, n_person, 1, n_person);
	free_dmatrix(chol, 1, n_person, 1, n_person);
}

/* omega calculates the variance-covariance matrix omega from
   the kinship/inbreeding coefficient matrix, the genotypes and the prespecified
   values of the genetic effec (gamma), additive variance (sigma_a) and error (sigma_e)
   HERE WE DO NOT CONSIDER GENETIC EFFECT - VALUE OF 0

   Our model is:
   MVN(effects of major genes, omega) = fixed function of  +  MVN(0, omega)
   the effect of the
   major genes

   gamma[marker][allele 1][allele 2] = effect
   mu[person]
   geno[person][locus][chrom]
   */
void omega_fun(double **omega, int n_person, double **kin, double add_var) {

	int i, j;

	for (i = 1; i <= n_person; i++) {
		for (j = 1; j <= n_person; j++) 
			omega[i][j] = add_var * kin[i][j];
	}
	
	// if (dom_var >= 10e-6) { // Can change that to ifelse above like (dom_var>=10e-6)? dom_var*delta_7[i][j]: 0
		// //printf("Dominant variance is included in phenotype model\n");
		// for (i = 1; i <= n_person; i++) {
			// for (j = 1; j <= n_person; j++) {
				// omega[i][j] += dom_var*delta_7[i][j];
			// }
		// }
	// }
}

void compute_trait_logistic(double *trait, double *pi, int n_person, long *seed) {
	int i;
	double tmp = 0;

	for (i = 1; i <= n_person; i++) {
		// Get pi and generate Y replicate
		tmp = exp(pi[i]) / (1.0 + exp(pi[i]));
		trait[i] = (rand1(seed) < tmp) ? 1 : 0;
	}
}

// void compute_trait_liability(int *trait, double *L, int n_person, struct SIM_PARAM param, long *seed) {
	// int i;
	// double var = .97 + param.add_var + param.err_var;
	// if (param.dom_var >= 10e-6) var += param.dom_var;

	// double T = 1.843 + sqrt(var)*0.52;

	// for (i = 1; i <= n_person; i++) trait[i] = (L[i] > T) ? 1 : 0;
// }

//int countallele(int *genotype){
//	int allele_count = 0;
//	if ((genotype[0] == 1) && (genotype[1] == 1))
//		allele_count = 2;
//	else if ((genotype[0] == 2) && (genotype[1] == 2))
//		allele_count = 0;
//	else
//		allele_count = 1;
//	return allele_count;
//}


