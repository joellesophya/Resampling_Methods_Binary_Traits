#include "simulation.h"

#include <R.h>
#include <Rinternals.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "nrutil.h"

#include "datacheck.h"
#include "read.h"
#include "trait_sim.h"
#include "sim_pheno_geno.h"
#include "assotest.h"

SEXP VT_analysis(SEXP Yvec, SEXP agevec, SEXP sexvec, SEXP Zvec, SEXP varval, SEXP betavec);

SEXP VT_analysis(SEXP Yvec, SEXP agevec, SEXP sexvec, SEXP Zvec, SEXP varval, SEXP betavec) {

	// Reading info from R
	double *Y, *age, *sex, *Zcov, est_add_var, *est_beta;
	Y = REAL(Yvec);
	age = REAL(agevec);
	sex = REAL(sexvec);
	Zcov = REAL(Zvec);
	est_add_var = REAL(varval)[0];
	est_beta = REAL(betavec); // intercept and covariates
	
	struct FAMILY *family;
	struct DATA_STRUCT data_struct;
	struct SIM_PARAM param;
	struct RV_INFO rv;
	struct STATS stats;

	int i, j, k;
	int n_fam, n_person, n_founder_total, n_cov, n_set, n_rep_cases, perm_accept;

	// Get the initial seed for all the rand functions
	// long seed_rand = 8734568;
	long seed_rand = (long)time(NULL);
	long seed = (-1)*seed_rand;
	srand(seed_rand);

	int n_perm = 5*1e4, i_perm;
	
	format_checks("pedigree", &data_struct);
	data_struct.n_cov = length(betavec);
	n_cov = data_struct.n_cov;
	n_fam = data_struct.n_fam;
	data_struct.tol = NUMTOL;
	
	family = (struct FAMILY *)malloc(sizeof(struct FAMILY)*(n_fam + 1));
	
	readsimparam("parameters", &param); /* Reads file that contains info about number of markers per rv set */
	data_struct.n_marker = param.n_marker;

	//printf("\nWARNING: User should make sure these counts are as expected:\n");
	//printf("n_fam = %d\nn_RV_marker/site = %d\nn_total = %d\nn_cov = %d\n", data_struct.n_fam, data_struct.n_marker, data_struct.n_total, data_struct.n_cov);

	pednum("pedigree", family, n_cov, data_struct.n_marker); /*Counts n.founders and n.nonfounders*/
	readped("pedigree", family); // Reads pedigree info
	readkininb("kinship", family); // Reads kinship coef
	n_founder_total = family[1].n_founder*n_fam;
	n_person = family[1].n_person;

	// Allocate space for matrices that will store RV info 
	rv.n_set = 150; // Number of RV sets being simulated (each set will correspond to a statistic)
	n_set = rv.n_set;
	rv.haplo = imatrix(1, 10000, 1, param.n_sites); // There are 10000 lines in haplo.txt
	rv.lines_haplo = imatrix(1, rv.n_set, 1, 2 * n_founder_total);
	rv.poly_rv = imatrix(1, n_set, 1, data_struct.n_marker);

	// Reads info about the haplotypes for founders and stores in integer matrix 'haplotype'
	readhaplo("../haplo.txt", rv.haplo, param.n_sites);

	// For OLS estimates used in VT
	data_struct.n_cases = ivector(1, n_perm + 1);
	
	// For VT statistic
	stats.MAF = dvector(1, data_struct.n_marker);
	stats.sorted_MAF = dvector(1, data_struct.n_marker);
	stats.unique_MAF = ivector(1, data_struct.n_marker);
	stats.threshold = dvector(1, data_struct.n_marker);
	stats.zVT = dvector(1, n_perm + 1);
	stats.mles = dvector(1, n_cov + 1); // Stores [beta_hat, add_var]
	stats.cumsum_num = dvector(1, n_perm + 1);
	stats.pvals = dvector(1, n_set);
	SEXP pvals_out = PROTECT(allocVector(REALSXP, n_set));

	// For trait matrix (holds Y and Y_pi for each family) -- indID is col
	for (i = 1; i <= n_fam; i++)	family[i].trait = dmatrix(1, n_perm + 1, 1, n_person);
	
	// Get trait, mle estimates and covariate
	for (i = 1, k = 0, data_struct.n_cases[1] = 0; i <= n_fam; i++){
		for (j = 1; j <= n_person; j++, k++){
			family[i].trait[1][j] = Y[k];
			family[i].cov[j][2] = age[k];
			family[i].cov[j][3] = sex[k];
			family[i].cov[j][n_cov] = Zcov[k];
			if (family[i].trait[1][j] == 1) data_struct.n_cases[1]++;
		}
	}
	for (i = 1; i <= n_cov; i++) stats.mles[i] = est_beta[i - 1];
	stats.mles[n_cov + 1] = est_add_var;
	
	// Simulate the replicates with conditioning on number of cases in sample
	i_perm = 1, i = 1;
	while(i_perm <= n_perm){
		for (n_rep_cases = 0, j = 1; j <= n_fam; j++) {
			run_trait_sim(family[j], stats.mles, n_cov, i_perm, &seed);	
			for (i = 1; i <= n_person; i++) n_rep_cases += (family[j].trait[i_perm + 1][i] == 1)? 1 : 0;
		}
		if(abs(n_rep_cases - data_struct.n_cases[1]) <= (0.05 * data_struct.n_total)){
			data_struct.n_cases[++i_perm] = n_rep_cases;
		}
		if (i++ > (n_perm / .01)){
			for (j = 0; j < n_set; j++) REAL(pvals_out)[j] = -1;
		}
	}
	
	if (i_perm > n_perm){  // Simulated enough Y replicates
		// Simulate genotype
		sim_pheno_geno(family, rv, data_struct, param.n_sites, &seed);
		
		// Get VT stat p-values for all RV sites
		VT(family, rv, n_set, data_struct, stats, n_perm, &seed);
		
		for (i = 1; i <= n_set; i++){
			est_add_var = stats.pvals[i];
			REAL(pvals_out)[i-1] = est_add_var;
		}
	}
	
	//	/////////////////////////////////////////////////////
	//	////////////////// Free memory //////////////////////
	//	/////////////////////////////////////////////////////
	for (i = 1; i <= n_fam; i++) {
		n_person = family[i].n_person;

		free_imatrix(family[i].parent, 1, n_person, 0, 1);
		free_dmatrix(family[i].phi, 1, n_person, 1, n_person);
		free_imatrix(family[i].genos, 1, n_person, 0, 1);
		free_dmatrix(family[i].cov, 1, n_person, 1, n_cov);
		free_dmatrix(family[i].trait, 1, n_perm + 1, 1, n_person);
	}

	free_ivector(data_struct.n_cases, 1, n_perm + 1);
	
	free_imatrix(rv.haplo, 1, 10000, 1, param.n_sites);
	free_imatrix(rv.lines_haplo, 1, rv.n_set, 1, 2 * n_founder_total);
	free_imatrix(rv.poly_rv, 1, rv.n_set, 1, data_struct.n_marker);

	free_dvector(stats.threshold, 1, data_struct.n_marker);
	free_dvector(stats.zVT, 1, n_perm + 1);
	free_dvector(stats.MAF, 1, data_struct.n_marker);
	free_dvector(stats.sorted_MAF, 1, data_struct.n_marker);
	free_ivector(stats.unique_MAF, 1, data_struct.n_marker);
	free_dvector(stats.mles, 1, n_cov + 1);
	free_dvector(stats.cumsum_num, 1, n_perm + 1);
	free_dvector(stats.pvals, 1, n_set);

	free(family);
	
	UNPROTECT(1);
	return pvals_out;
} /* end main */
