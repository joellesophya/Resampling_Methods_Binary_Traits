

// [[1]]
     // alpha Err_rate SE     p_value
// [1,] 0.005 0.005    0.0007 1      
// [2,] 0.01  0.0099   0.001  0.9599 
// [3,] 0.05  0.0519   0.0022 0.396  

// [[2]]
     // alpha Err_rate SE     p_value
// [1,] 0.005 0.0051   0.0007 0.9435 
// [2,] 0.01  0.0103   0.001  0.8016 
// [3,] 0.05  0.0518   0.0022 0.422  

// [[3]]
     // alpha Err_rate SE     p_value
// [1,] 0.005 0.0049   0.0007 0.9435 
// [2,] 0.01  0.0099   0.001  0.9599 
// [3,] 0.05  0.0516   0.0022 0.477 

#include "simulation.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "nrutil.h"

#include "sim_pheno_geno.h"

/* main */
int main(int argc, char **argv) {

	struct FAMILY *family;
	struct DATA_STRUCT data_struct;
	struct SIM_PARAM param;
	struct RV_INFO rv;
	struct STATS stats;
	struct PHI_SVD phi_svd;
	struct TMP_PERM tmp_perm;

	double alpha, p_val, tmp;
	double *p_vals; // Vector that stores p_values
	char perm_filename[20];

	int i, j, k;
	int n_fam, n_person, n_founder_total, n_cov, n_set;

	// For VT
	double *MAF, *sorted_MAF;
	int *unique_MAF;
	double *threshold, *z;
	double *VT_obs;

	// Get the initial seed for all the rand functions
	//long seed_rand = 8734568;
	long seed_rand = (long)time(NULL);
	long seed = (-1)*seed_rand;
	srand(seed_rand);

	// Define number of phenotype relaizations and number of permutation replicates
	int n_rep = 100, i_rep;
	int perm_num[3] = { 15000, 20000, 25000 };
	//int perm_num[3] = { 10, 20, 30 };
	int n_perm = perm_num[2], i_perm, iter_perm = 0;
	int gen_model = 1; // Which model to use for generating phenotype? 1=logistic (default) vs. 2=liability threshold

	// File to store p_value info
	FILE *p_vals_out;
	// Store the seed info in a file
	FILE *seed_file;
	seed_file = fopen("seed.txt", "w");
	fprintf(seed_file, "%ld\n", seed_rand);
	fclose(seed_file);

	//** For each Y realization **//
	for (i_rep = 1; i_rep <= n_rep; i_rep++){

		if (i_rep == 1){
			format_checks("pedigree", &data_struct);
			n_fam = data_struct.n_fam;
			data_struct.tol = NUMTOL;
			n_cov = data_struct.n_cov;

			family = (struct FAMILY *)malloc(sizeof(struct FAMILY)*(n_fam + 1));

			readsimparam("parameters", &param); /* Reads file that contains info about number of parameters */
			data_struct.n_marker = param.n_marker;

			printf("\nWARNING: User should make sure these counts are as expected:\n");
			printf("n_fam = %d\nn_RV_marker/site = %d\nn_total = %d\nn_cov = %d\n", data_struct.n_fam, data_struct.n_marker, data_struct.n_total, data_struct.n_cov);

			pednum("pedigree", family, n_cov, data_struct.n_marker); /*Counts n.founders and n.nonfounders*/
			readped("pedigree", family); /*reads pedigree info*/
			readkininb("kinship", family); /*Reads kinship & inbreeding coef*/
			n_founder_total = family[1].n_founder*n_fam;

			// Allocate space for matrices that will store RV info 
			rv.n_set = 100; // Number of RV sets being simulated (each set will correspond to a statistic)
			n_set = rv.n_set;
			rv.haplo = imatrix(1, 10000, 1, param.n_sites); // There are 10000 lines in haplo.txt
			rv.lines_haplo = imatrix(1, rv.n_set, 1, 2 * n_founder_total);
			rv.poly_rv = imatrix(1, n_set, 1, data_struct.n_marker);

			// Reads info about the haplotypes for founders and stores in integer matrix 'haplotype'
			readhaplo("haplo.txt", rv.haplo, param.n_sites);

			// For VT statistics observed
			VT_obs = dvector(1, n_set);
			MAF = dvector(1, data_struct.n_marker);
			sorted_MAF = dvector(0, data_struct.n_marker - 1);
			unique_MAF = ivector(1, data_struct.n_marker);
			threshold = malloc((data_struct.n_marker + 1)*sizeof(double));
			z = malloc((data_struct.n_marker + 1)*sizeof(double));

			phi_svd.u = dmatrix(1, family[1].n_person, 1, family[1].n_person); // Assume all families have the same number of individuals (i.e. n_person)
			phi_svd.lambda = dvector(1, family[1].n_person);
			phi_svd.sigma = dmatrix(1, family[1].n_person, 1, family[1].n_person);
			compute_phi_svd(family[1], phi_svd);

			tmp_perm.delta = dvector(1, data_struct.n_total - data_struct.n_cov);
			tmp_perm.delta_perm = dvector(1, data_struct.n_total - data_struct.n_cov);
			tmp_perm.v1 = dmatrix(1, data_struct.n_total, 1, data_struct.n_total - data_struct.n_cov);
			tmp_perm.cholesky = dmatrix(1, family[1].n_person, 1, family[1].n_person);
			tmp_perm.v1Delta = dvector(1, data_struct.n_total);// Used in the permute function
			tmp_perm.v1Delta_f = dmatrix(1, n_fam, 1, family[1].n_person);

			p_vals = dvector(1, n_set); // Vector that stores p_values
		}

		// Simulate phenotype and genotype
		sim_pheno_geno(family, rv, data_struct, param, gen_model, &seed);

		// Get VT stat for each set of RV sites for Y realization
		for (i = 1; i <= n_set; i++){
			VT(family, rv, i, data_struct, &stats, MAF, sorted_MAF, unique_MAF, threshold, z, 0);
			VT_obs[i] = fabs(stats.vt);
		}

		////////////////////////////////////////////////////////
		/*--------Marker analysis and getting null VCs--------*/
		////////////////////////////////////////////////////////
		compute_gee_res(family, phi_svd, data_struct);

		////////////////////////////////////////////////////////////////
		///                 PERMUTATION BASED METHOD                  //
		////////////////////////////////////////////////////////////////
		// Generate delta vector that is common to all permutation replicates (before shuffling entries of delta -- see notes)
		permute_pre(family, data_struct, phi_svd, tmp_perm);

		// For p-values
		for (i = 1; i <= n_set; i++) p_vals[i] = 0;

		iter_perm = 0;
		for (i_perm = 1; i_perm <= n_perm; i_perm++){
			// Get permutation replicate and store it in family.y_perm
			permute(family, data_struct, tmp_perm);

			// Get VT stats for each rv set
			for (i = 1; i <= n_set; i++){
				VT(family, rv, i, data_struct, &stats, MAF, sorted_MAF, unique_MAF, threshold, z, 1);
				if (fabs(stats.vt) >= VT_obs[i]) p_vals[i] += 1;
			}

			if (i_perm == perm_num[iter_perm]){
				sprintf(perm_filename, "pvals_%d.txt", i_perm);
				p_vals_out = fopen(perm_filename, (i_rep == 1) ? "w" : "a");

				for (i = 1; i < n_set; i++) fprintf(p_vals_out, "%.7f ", (p_vals[i] + 1.0) / (i_perm + 1.0));
				fprintf(p_vals_out, "%.7f\n", (p_vals[n_set] + 1.0) / (i_perm + 1.0));

				fclose(p_vals_out);
				iter_perm++;
			}
		}
	}
	//	/////////////////////////////////////////////////////
	//	////////////////// Free memory //////////////////////
	//	/////////////////////////////////////////////////////
	for (i = 1; i <= n_fam; i++) {
		n_person = family[i].n_person;

		free_dmatrix(family[i].phi, 1, n_person, 1, n_person);
		free_dmatrix(family[i].delta_7, 1, n_person, 1, n_person);
		free_imatrix(family[i].genos, 1, n_person, 0, 1);
		free_dmatrix(family[i].cov, 1, n_person, 1, n_cov);
		free_imatrix(family[i].parent, 1, n_person, 0, 1);
		free_dvector(family[i].trait, 1, n_person);
		free_dvector(family[i].mu_hat, 1, n_person);
		free_dvector(family[i].y_perm, 1, n_person);
	}

	free_imatrix(rv.haplo, 1, 10000, 1, param.n_sites);
	free_imatrix(rv.lines_haplo, 1, rv.n_set, 1, 2 * n_founder_total);
	free_imatrix(rv.poly_rv, 1, rv.n_set, 1, data_struct.n_marker);

	free_dvector(VT_obs, 1, n_set);
	free(threshold);
	free(z);
	free_dvector(MAF, 1, data_struct.n_marker);
	free_dvector(sorted_MAF, 0, data_struct.n_marker - 1);
	free_ivector(unique_MAF, 1, data_struct.n_marker);

	free_dmatrix(phi_svd.u, 1, family[1].n_person, 1, family[1].n_person); // Assume all families have the same number of individuals (i.e. n_pheno)
	free_dvector(phi_svd.lambda, 1, family[1].n_person);
	free_dmatrix(phi_svd.sigma, 1, family[1].n_person, 1, family[1].n_person);

	free_dvector(tmp_perm.delta, 1, data_struct.n_total - data_struct.n_cov);
	free_dvector(tmp_perm.delta_perm, 1, data_struct.n_total - data_struct.n_cov);
	free_dmatrix(tmp_perm.v1, 1, data_struct.n_total, 1, data_struct.n_total - data_struct.n_cov);
	free_dmatrix(tmp_perm.cholesky, 1, family[1].n_person, 1, family[1].n_person);
	free_dvector(tmp_perm.v1Delta, 1, data_struct.n_total);// Used in the permute function
	free_dmatrix(tmp_perm.v1Delta_f, 1, n_fam, 1, family[1].n_person);

	free(family);

	free_dvector(p_vals, 1, n_set); // Vector that stores p_values

	return 0;
} /* end main */