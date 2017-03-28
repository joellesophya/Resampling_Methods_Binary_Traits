#include "simulation.h"

#include <R.h>
#include <Rdefines.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "nrutil.h"

#include "datacheck.h"
#include "read.h"
#include "permute_sim.h"
#include "sim_pheno_geno.h"
#include "assotest.h"

SEXP VT_analysis(SEXP Yvec, SEXP ncov, SEXP mu_Y, SEXP delt, SEXP GCVm, SEXP nperm, SEXP method_name, SEXP Yperm_type, SEXP nRVs, SEXP nsets, SEXP n_markers);

SEXP VT_analysis(SEXP Yvec, SEXP ncov, SEXP mu_Y, SEXP delt, SEXP GCVm, SEXP nperm, SEXP method_name, SEXP Yperm_type, SEXP nRVs, SEXP nsets, SEXP n_markers) {
	// Reading info from R
	int n_cov, n_perm, convert_to_bin, n_sites, n_set;
	n_cov = INTEGER(ncov)[0]; // # of covariates in EE method
	n_perm = INTEGER(nperm)[0]; // number of trait replicates to generate 
	convert_to_bin = INTEGER(Yperm_type)[0]; // should trait replicate be converted to binary?
	
	// Determine if using VT or WS test statistic
	PROTECT(method_name = AS_CHARACTER(method_name));
	char *method[length(method_name)];
	method[0] = R_alloc(strlen(CHAR(STRING_ELT(method_name, 0))), sizeof(char)); 
	strcpy(method[0], CHAR(STRING_ELT(method_name, 0)));
	if( strcmp(method[0], "VT") != 0 ) convert_to_bin = 1;   // Convert replicates to binary if using WST
	
	struct FAMILY *family;
	struct DATA_STRUCT data_struct;
	struct RV_INFO rv;
	struct STATS stats;

	int i, j, fam, index_tot;
	int n_total, n_fam, n_person, n_founder_total;

	// Get the initial seed for all the rand functions
	// long seed_rand = 8734568;
	long seed_rand = (long)time(NULL);
	long seed = (-1)*seed_rand;
	srand(seed_rand);

	format_checks("pedigree", &data_struct);
	data_struct.n_cov = n_cov;
	n_fam = data_struct.n_fam;
	n_total = data_struct.n_total;
	data_struct.tol = NUMTOL;
	// printf(" There are %d individuals among %d families and %d covariates\n", n_total, n_fam, n_cov);
	
	family = (struct FAMILY *)malloc(sizeof(struct FAMILY)*(n_fam + 1));
	
	n_sites = INTEGER(nRVs)[0];
	n_set = INTEGER(nsets)[0]; // Number of RV sets being simulated (each set will correspond to a statistic)
	data_struct.n_marker = INTEGER(n_markers)[0];
	pednum("pedigree", family); /*Counts n.founders and n.nonfounders*/
	readped("pedigree", family); // Reads pedigree info
	n_founder_total = family[1].n_founder * n_fam; //same pedigree structure
	n_person = family[1].n_person;

	// Allocate space for matrices that will store RV info 
	rv.n_set = n_set; 
	rv.haplo = imatrix(1, 10000, 1, n_sites); // There are 10000 lines in haplo.txt
	rv.lines_haplo = imatrix(1, n_set, 1, 2 * n_founder_total);
	rv.poly_rv = imatrix(1, n_set, 1, data_struct.n_marker);

	// Reads info about the haplotypes for founders and stores in integer matrix 'haplotype'
	readhaplo("../haplo.txt", rv.haplo, n_sites);

	// For VT statistic
	data_struct.sumsY = dvector(1, 2 * (n_perm + 1)); // Stores {sum(Y), sum(Y^2)}
	stats.delta = dvector(1, n_total - n_cov); 
	stats.GCV_mat = dmatrix(1, n_total, 1, n_total - n_cov); 
	stats.pvals = dvector(1, n_set);
	SEXP pvals_out = PROTECT(allocVector(REALSXP, n_set));

	// For trait matrix (holds Y and Y_pi for each family) -- indID is col
	for (i = 1; i <= n_fam; i++) family[i].trait = dmatrix(1, n_perm + 1, 1, n_person);
	
	// Get trait, covariate, pi_hat, delta vector & GCV matrix for EE-perm. method
	for (data_struct.sumsY[1] = data_struct.sumsY[2 + n_perm] = 0, fam = 1; fam <= n_fam; fam++){
		for (family[fam].naff = 0, i = 1; i <= n_person; i++){
			index_tot = (fam - 1) * n_person + i; // Index in whole sample (from 1 to n_tot)
			// Store trait values and mu_Y_hat
			family[fam].trait[1][i] = REAL(Yvec)[index_tot - 1]; // Y[0] to Y[n_tot-1])
			family[fam].naff += family[fam].trait[1][i];
			data_struct.sumsY[1] += REAL(Yvec)[index_tot - 1];
			data_struct.sumsY[2 + n_perm] += REAL(Yvec)[index_tot - 1] * REAL(Yvec)[index_tot - 1];
			family[fam].mu_hat[i] = REAL(mu_Y)[index_tot - 1];
			for(j = 1; j <= (n_total - n_cov); j++){
				// Store GCV matrix -- read by row
				stats.GCV_mat[index_tot][j] = REAL(GCVm)[(index_tot - 1) * (n_total - n_cov) + j - 1];
				// Store delta vector
				if((fam == 1) && (i == 1)) stats.delta[j] = REAL(delt)[j - 1];
			}
		}
	}
	
	// Simulate genotype
	sim_pheno_geno(family, n_fam, rv, data_struct.n_marker, n_sites, &seed);
	// Simulate the replicates
	data_struct.sumsY = permute(family, data_struct, stats, n_perm, convert_to_bin);
	
	// Get p-values for all RV sites
	if( strcmp(method[0], "VT") == 0 )	VT(family, rv, n_set, data_struct, stats.pvals, n_perm, &seed);
	else ws(family, rv, n_set, data_struct, stats.pvals, n_perm, &seed);
	
	for (i = 1; i <= n_set; i++) REAL(pvals_out)[i-1] = stats.pvals[i];
	
	//	/////////////////////////////////////////////////////
	//	////////////////// Free memory //////////////////////
	//	/////////////////////////////////////////////////////
	for (i = 1; i <= n_fam; i++) {
		n_person = family[i].n_person;

		free_imatrix(family[i].parent, 1, n_person, 0, 1);
		// free_dmatrix(family[i].phi, 1, n_person, 1, n_person);
		free_imatrix(family[i].genos, 1, n_person, 0, 1);
		free_dmatrix(family[i].trait, 1, n_perm + 1, 1, n_person);
		free_dvector(family[i].mu_hat, 1, n_person);
	}

	free_dvector(data_struct.sumsY, 1, 2 * (n_perm + 1));
	
	free_imatrix(rv.haplo, 1, 10000, 1, n_sites);
	free_imatrix(rv.lines_haplo, 1, n_set, 1, 2 * n_founder_total);
	free_imatrix(rv.poly_rv, 1, n_set, 1, data_struct.n_marker);

	free_dvector(stats.delta, 1, n_total - n_cov);
	free_dmatrix(stats.GCV_mat, 1, n_total, 1, n_total - n_cov); 
	free_dvector(stats.pvals, 1, n_set);

	free(family);
	
	UNPROTECT(2);
	return pvals_out;
} /* end main */
