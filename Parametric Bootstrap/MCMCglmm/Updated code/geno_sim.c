#include <stdio.h>
#include <stdlib.h>
#include "sim_pheno_geno.h"
#include "geno_sim.h"
#include "nrutil.h"
#include "rand1.h"

void drop(int nfou, int nnon, int **parentmat, int **geno, long *seed) {
	int i, parent, trchrom;
	int npeople, fcount;

	npeople = nfou + nnon;
	for (fcount = 0, i = 1; i <= npeople; i++) {
		/* If the person is a founder, enter founder genotypes in geno array */
		if (parentmat[i][0] == 0 && parentmat[i][1] == 0) {
			fcount++;
			// Assign number 1 to (2*n_founder) to each haplotype of the founders
			geno[i][0] = 2 * fcount - 1;
			geno[i][1] = 2 * fcount;

		} /* else if both parents are in the pedigree drop alleles down the pedigree */
		else if (parentmat[i][0] != 0 && parentmat[i][1] != 0) {
			for (parent = 0; parent < 2; parent++){
				/* Randomly choose chrom. 0 or 1 to transmit to kid at each locus (Assume markers are on different chromosomes*/
				trchrom = (int)(rand1(seed) * 2);
				// For COSI haplotype
				geno[i][parent] = geno[parentmat[i][parent]][trchrom];

			}
		} /* else if only one parent is in the pedigree, write out an error message*/
		else {
			fprintf(stderr, "%dth person has only one unknown parent. Illegal.\n", i);
			exit(2);
		}
	}
}

void codealleles(int n_founder_allele, int *state, long *seed) {
	int i, assigned;
	double frac, MAF = 0.3; // Common SNPs with MAF of .3

	for (i = 1; i <= n_founder_allele; i++){
		assigned = 0;
		frac = rand1(seed);
		if (frac < MAF){
			state[i] = 1; // 1 is the minor allele
			assigned = 1;
		}
		else{
			state[i] = 2; // 2 is the major allele
			assigned = 1;
		}
		if (assigned == 0) {
			printf("not assigned\n");
		}
	}
}

void recode(int n_person, int *state, int **geno_in, int **geno_out) {
	int i;

	for (i = 1; i <= n_person; i++) {
		geno_out[i][0] = state[geno_in[i][0]];
		geno_out[i][1] = state[geno_in[i][1]];
	}
}

void run_geno_sim(struct FAMILY family, int n_marker, long *seed) {
	int n_founder = family.n_founder;
	int n_non_founder = family.n_non_founder;

	int **parent = family.parent;
	int **geno_tmp = family.genos;

	drop(n_founder, n_non_founder, parent, geno_tmp, seed);

}
