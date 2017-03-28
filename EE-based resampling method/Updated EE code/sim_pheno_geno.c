#include "sim_pheno_geno.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "nrutil.h"
#include "geno_sim.h"

//////////////////////////////////////////////////
//////// Simulate phenotype and genotype /////////
//////////////////////////////////////////////////

void sim_pheno_geno(struct FAMILY *family, int n_fam, struct RV_INFO rv, int n_marker, int n_sites, long *seed) {
	int fam, n_founder_total = family[1].n_founder * n_fam; // Assume same number of founders accross families

	// Choose founder haplotypes and polymorphic sites
	store_haploline(rv, n_founder_total, n_sites, n_marker);

	// Generate genotype data = Assign all individuals two founder haplotypes
	for (fam = 1; fam <= n_fam; fam++) run_geno_sim(family[fam], n_marker, seed); 
}

// Store haplotype line information for founders' haplotypes in pedigree and also info about first 50 polymorphic sites
void store_haploline(struct RV_INFO rv, int n_founders, int n_sites, int n_marker){
	int i, k, i_set, n_set = rv.n_set;
	int lines_in_haplo = 10000;
	int **lines = rv.lines_haplo;
	int **haplo = rv.haplo;
	int **poly_rvs = rv.poly_rv;

	// Choose (2*n_founder) numbers without replacement from 1-10000 for each set of rv's
	// Those will be the haplotypes for the founders in sample
	int *lines_haplo = ivector(1, lines_in_haplo);
	int *scramble_cols_haplo = ivector(1, n_sites);

	for (i = 1; i <= lines_in_haplo; i++){
		lines_haplo[i] = i; // Vector [1 2 ... 9999 10000]
		if(i <= n_sites) scramble_cols_haplo[i] = i;
	}
	for (i_set = 1; i_set <= n_set; i_set++){
		shuffle_int(lines_haplo, lines_in_haplo);
		for (i = 1; i <= (2 * n_founders); i++)	lines[i_set][i] = lines_haplo[i];
	}

	// Choose first n_marker polymorphic sites 
	// For each rv set, loop though all RV sites (= col in haplo matrix)
	int n_poly, col_haplo; // Keep count of identity and number of polymorphic sites found
	int allele_haplo1; // used to check polymorphicity

	for (i_set = 1; i_set <= n_set; i_set++){
		n_poly = 1; // Counts # of polymorphic sites
		shuffle_int(scramble_cols_haplo, n_sites);
		for (i = 1; i <= n_sites; i++){
			// For each site, check if it is polymorphic accross the (2*n_founders) haplotypes
			col_haplo = scramble_cols_haplo[i];
			allele_haplo1 = haplo[lines[i_set][1]][col_haplo];
			for (k = 2; k <= (2 * n_founders); k++){
				if (allele_haplo1 != haplo[lines[i_set][k]][col_haplo]){
					poly_rvs[i_set][n_poly++] = col_haplo; // Store RV#
					break;
				}
			}
			// If enough polymorphic markers have been found, break the initial for loop
			if (n_poly > n_marker)
				break;
			// If went through all sites and did not find enough polymorphic sites
			if ((i == n_sites) && (n_poly <= n_marker)){
				printf("Not enough polymorphic sites in founder haplotypes selected\n");
				exit(1);
			}
		}
	}
	free_ivector(lines_haplo, 1, lines_in_haplo);
	free_ivector(scramble_cols_haplo, 1, n_sites);
}

static int rand_int(int n) {
	int limit = RAND_MAX - RAND_MAX % n;
	int rnd;
	do {
		rnd = rand();
	} while (rnd >= limit);
	return rnd % n; // Value between 0 and (n-1)
}

// n is length of the array
void shuffle(double *array, int n) {
	int i, j;
	double tmp;
	for (i = n; i > 1; i--) {
		j = rand_int(i) + 1;
		tmp = array[j];
		array[j] = array[i];
		array[i] = tmp;
	}
}

// n is length of the array of integers
void shuffle_int(int *array, int n) {
	int i, j, tmp;
	for (i = n; i > 1; i--) {
		j = rand_int(i) + 1;
		tmp = array[j];
		array[j] = array[i];
		array[i] = tmp;
	}
}
