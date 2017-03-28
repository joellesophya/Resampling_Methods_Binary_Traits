// # 100 Y rep, 15000 Y perm, sigma_a^2=4 using EE-logistic model

// > get_type1err(p); get_type1err(pASTOR)
     // alpha Err_rate SE     p_value
// [1,] 0.005 0.0047   0.0003 0.428  
// [2,] 0.01  0.0094   0.0004 0.2164 
// [3,] 0.05  0.0463   0.0009 0.0001 
     // alpha Err_rate SE     p_value
// [1,] 0.005 0.0053   0.0003 0.428  
// [2,] 0.01  0.0103   0.0005 0.544  
// [3,] 0.05  0.0498   0.001  0.8777 
// Estimating equation method
#include "simulation.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "nrutil.h"

#include "sim_pheno_geno.h"
#include "mle.h"
struct OPT_STRUCT f_opt_struct;
void init_opt_struct(struct FAMILY *family, struct FAMILY_SVD *svd, struct DATA_STRUCT data_struct) {
	f_opt_struct.family = family;
	f_opt_struct.svd = svd;
	f_opt_struct.data_struct = data_struct;
}

/* main */
int main(int argc, char **argv) {

	struct FAMILY *family;
	struct DATA_STRUCT data_struct;
	struct SIM_PARAM param;
	struct ASTOR astor;
	struct MLE_PARAM mle_initial;
	struct FAMILY_SVD *svd;
	struct PHI_SVD phi_svd;
	struct MLE_RESULTS mle_res;
	struct TMP_PERM tmp_perm;

	double alpha, p_val, tmp;
	double *v1Delta; // For obtaining the permuted phenotype vector
	double **v1Delta_f;
	double *p_vals; // Vector that stores p_values
	int i, j, k, n_fam, n_person, n_cov, nullmle = 0;

	/**************************/
	/********** Done *********/
	// Set the seed fo rand() to the time the app is ran
	long seed_rand = (long)time(NULL);
	long seed = (-1)*seed_rand;
	srand(seed_rand);
	printf("The initial seed used for Y replicate is: %ld\n", seed_rand);

	// Define number of phenotype relaizations and replicates
	int n_rep = 100, i_rep;
	int n_perm = 15000, i_perm;
	int gen_model = 1; // 1= Default is logistic (2 is liability threshold)

	// File to store p_value info
	FILE *ASTOR_stats;
	FILE *p_vals_out;

	//** For each Y realization **//
	for (i_rep = 1; i_rep <= n_rep; i_rep++){

		if (i_rep == 1){
			format_checks("pedigree", &data_struct);
			n_fam = data_struct.n_fam;
			data_struct.tol = NUMTOL;
			n_cov = data_struct.n_cov;

			family = (struct FAMILY *)malloc(sizeof(struct FAMILY)*(n_fam + 1));
			svd = (struct FAMILY_SVD *)malloc(sizeof(struct FAMILY_SVD)*(n_fam + 1));

			readsimparam("parameters", &param); /* Reads file that contains info about number of parameters */
			data_struct.n_marker = param.n_marker;

			printf("\nWARNING: User should make sure these counts are as expected:\n");
			printf("n_fam = %d\nn_marker = %d\nn_total = %d\nn_cov = %d\n", data_struct.n_fam, data_struct.n_marker, data_struct.n_total, data_struct.n_cov);

			pednum("pedigree", family, n_cov, data_struct.n_marker); /*Counts n.founders and n.nonfounders*/
			readped("pedigree", family); /*reads pedigree info*/
			readkininb("kinship", family); /*Reads kinship & inbreeding coef*/

			astor.C = dmatrix(1, family[1].n_person, 1, family[1].n_person);
			astor.sigma_0 = dvector(1, data_struct.n_marker);
			astor.stats = dvector(1, data_struct.n_marker);
			astor.perm_stats = dvector(1, data_struct.n_marker);

			phi_svd.u = dmatrix(1, family[1].n_person, 1, family[1].n_person); // Assume all families have the same number of individuals (i.e. n_pheno)
			phi_svd.lambda = dvector(1, family[1].n_person);
			phi_svd.sigma = dmatrix(1, family[1].n_person, 1, family[1].n_person);

			/*for (i = 1; i <= n_fam; i++){
				svd[i].trait_u = dvector(1, family[i].n_person);
				svd[i].cov_u = dmatrix(1, family[i].n_person, 1, data_struct.n_cov);
				svd[i].lambda = dvector(1, family[i].n_person);
				}*/

			tmp_perm.delta = dvector(1, data_struct.n_total - data_struct.n_cov);
			tmp_perm.v1 = dmatrix(1, data_struct.n_total, 1, data_struct.n_total - data_struct.n_cov);
			tmp_perm.cholesky = dmatrix(1, family[1].n_person, 1, family[1].n_person);
			for (i = 1; i <= n_person; i++){
				phi_svd.lambda[i] = 0;
				for (j = 1; j <= n_person; j++) {
					phi_svd.u[i][j] = 0;
					phi_svd.sigma[i][j] = 0;
					tmp_perm.cholesky[i][j] = 0;
				}
			}
			compute_phi_svd(family[1], phi_svd);

			v1Delta = dvector(1, data_struct.n_total);// Used in the permute function
			v1Delta_f = dmatrix(1, n_fam, 1, family[1].n_person);

			p_vals = dvector(1, data_struct.n_marker); // Vector that stores p_values
		}

		// Simulate phenotype and genotype
		sim_pheno_geno(family, data_struct, param, gen_model, &seed);

		////////////////////////////////////////////////////////
		/*------------------ASTOR ANALYSIS---------------------*/
		////////////////////////////////////////////////////////
		if (i_rep % 20 == 0) printf("\nStarting analysis for rep %d: ", i_rep);

		// Getting the MLE for beta, (not sigma_a^2 and sigma_e^2)
		mle_initial.tol = NUMTOL;

		mle_res.var_beta_hat = dvector(1, data_struct.n_cov);
		mle_res.beta_hat = dvector(1, data_struct.n_cov);
		mle_res.beta_ols = dvector(1, data_struct.n_cov);

		estimate_mle_ols(&mle_res, family, data_struct, 2, 0);

		// Estimate the nuisance parameters under the null hypothesis
		compute_gee_res(family, phi_svd, data_struct);

		compute_ASTOR(family, data_struct, mle_res, astor);
		///////////////////////////////////////////////////////////////
		//                 PERMUTATION BASED METHOD                  //
		///////////////////////////////////////////////////////////////
		permute_pre(family, data_struct, phi_svd, tmp_perm);
		for (i = 1; i <= data_struct.n_marker; i++) p_vals[i] = 0;

		for (i_perm = 1; i_perm <= n_perm; i_perm++){
			//if (i_perm % 1000 == 0) printf("%d ", i_perm);
			// Get permutation replicate and store it in family.y_perm
			permute(family, data_struct, tmp_perm, v1Delta, v1Delta_f);

			compute_ASTOR_perm(family, data_struct, mle_res, astor);
			for (i = 1; i <= data_struct.n_marker; i++){
				/*for (j = 1; j <= data_struct.n_marker; j++){
					if (astor.perm_stats[j] > astor.stats[i]) p_vals[i] += 1;
					else if (astor.perm_stats[j] == astor.stats[i]) p_vals[i] += ((rand1(&seed) < .5) ? 1 : 0);
					}*/

				if (astor.perm_stats[i] > astor.stats[i]) p_vals[i] += 1;
				else if (astor.perm_stats[i] == astor.stats[i]) p_vals[i] += ((rand1(&seed) < .5) ? 1 : 0);

			}
		} //printf("\n");

		//	Store ASTOR statistics as well as permutation based p-values
		ASTOR_stats = fopen("ASTOR_statistics.txt", (i_rep == 1) ? "w" : "a");
		p_vals_out = fopen("sim_pvals.txt", (i_rep == 1) ? "w" : "a");

		for (i = 1; i < data_struct.n_marker; i++){
			//Compute p-value and store in in sim_pvals file
			p_vals[i] = (p_vals[i] + 1) / (n_perm + 1); //(n_perm*data_struct.n_marker + 1);
			fprintf(ASTOR_stats, "%f ", astor.stats[i]);
			fprintf(p_vals_out, "%f ", p_vals[i]);
		}
		fprintf(ASTOR_stats, "%f\n", astor.stats[data_struct.n_marker]);
		p_vals[data_struct.n_marker] = (p_vals[data_struct.n_marker] + 1) / (n_perm + 1); //(n_perm*data_struct.n_marker + 1);
		fprintf(p_vals_out, "%f\n", p_vals[data_struct.n_marker]);
		fclose(ASTOR_stats);
		fclose(p_vals_out);

		// Free memory for mle vectors (perhaps change this to do once within whole simulation???)
		free_dvector(mle_res.var_beta_hat, 1, data_struct.n_cov);
		free_dvector(mle_res.beta_hat, 1, data_struct.n_cov);
		free_dvector(mle_res.beta_ols, 1, data_struct.n_cov);
	}
	//	/////////////////////////////////////////////////////
	//	////////////////// Free memory //////////////////////
	//	/////////////////////////////////////////////////////
	for (i = 1; i <= n_fam; i++) {
		n_person = family[i].n_person;
		/*if (svd[i].n_ind != 0) {
			free_dvector(svd[i].trait_u, 1, n_person);
			free_dmatrix(svd[i].cov_u, 1, n_person, 1, n_cov);
			free_dvector(svd[i].lambda, 1, n_person);
			}*/
		free_dmatrix(family[i].phi, 1, n_person, 1, n_person);
		free_dmatrix(family[i].delta_7, 1, n_person, 1, n_person);
		free_i3tensor(family[i].genos, 1, data_struct.n_marker, 1, n_person, 0, 1);
		free_dmatrix(family[i].cov, 1, n_person, 1, n_cov);
		free_imatrix(family[i].parent, 1, n_person, 0, 1);
		free_dvector(family[i].trait, 1, n_person);
		free_dvector(family[i].y_perm, 1, n_person);
		free_dvector(family[i].v1_d, 1, n_person);
		free_dvector(family[i].mu_hat, 1, n_person);
		free_dmatrix(family[i].omega, 1, n_person, 1, n_person);
	}
	free_dmatrix(astor.C, 1, family[1].n_person, 1, family[1].n_person);
	free_dvector(astor.sigma_0, 1, data_struct.n_marker);
	free_dvector(astor.stats, 1, data_struct.n_marker);
	free_dvector(astor.perm_stats, 1, data_struct.n_marker);

	free_dmatrix(phi_svd.u, 1, family[1].n_person, 1, family[1].n_person); // Assume all families have the same number of individuals (i.e. n_pheno)
	free_dvector(phi_svd.lambda, 1, family[1].n_person);
	free_dmatrix(phi_svd.sigma, 1, family[1].n_person, 1, family[1].n_person);

	free_dvector(tmp_perm.delta, 1, data_struct.n_total - data_struct.n_cov);
	free_dmatrix(tmp_perm.v1, 1, data_struct.n_total, 1, data_struct.n_total - data_struct.n_cov);
	free_dmatrix(tmp_perm.cholesky, 1, family[1].n_person, 1, family[1].n_person);

	free_dvector(v1Delta, 1, data_struct.n_total);
	free_dmatrix(v1Delta_f, 1, n_fam, 1, family[1].n_person);

	free_dvector(p_vals, 1, data_struct.n_marker); // Vector that stores p_values

	free(family);
	free(svd);

	return 0;
} /* end main */