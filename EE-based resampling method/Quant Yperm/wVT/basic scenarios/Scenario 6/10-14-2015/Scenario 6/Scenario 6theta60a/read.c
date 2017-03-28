#include "read.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "nrutil.h"


// Read in parameter file
void readsimparam(char *filename, struct SIM_PARAM *param) {

	FILE *paramfile;
	paramfile = fopen(filename, "r");

	if (paramfile == NULL) {
		printf("Can't open parameters file\n");
		exit(1);
	}

	fscanf(paramfile, "%d ", &param->n_sites);
	fscanf(paramfile, "%d ", &param->n_marker);
	fscanf(paramfile, "%lf ", &param->add_var);
	fscanf(paramfile, "%lf ", &param->err_var);
	fscanf(paramfile, "%lf ", &param->mean); // This is the intercept term
	fscanf(paramfile, "%lf ", &param->dom_var);
	fscanf(paramfile, "%lf ", &param->sex_effect);
	fscanf(paramfile, "%lf ", &param->age_effect);
	fscanf(paramfile, "%d ", &param->cov_type);
	fscanf(paramfile, "%d ", &param->err_type);
	fscanf(paramfile, "%d ", &param->fam0);

	fclose(paramfile);
}

/* This function (pednum) counts the number of founders and non-founders in a pedigree.
   It takes in pointers to the variables: n_founder and n_non_founder and then uses
   temp variables to go through the whole file and count lines of specific format, that
   either have both parents zero (i.e. founder) or both parents non-zero (i.e. non-founder).
   It writes out an error message if one parent is zero while the other is non-zero.
   */

void pednum(char* filename, struct FAMILY *family, int n_cov, int n_marker) {

	int n_founder = 0;
	int n_non_founder = 0;
	int n_person = 0;

	/* Start with family ID=1*/
	int fam, fam0 = 1, ind, fa, mo, sex, age;

	FILE *pedfile;

	pedfile = fopen(filename, "r");

	if (pedfile == NULL) {
		printf("Can't open pedigree file\n");
		exit(1);
	}

	while (!feof(pedfile)) {
		fscanf(pedfile, "%d %d %d %d %d %d \n", &fam, &ind, &fa, &mo, &sex, &age); /*Changed it so could correctly read the pedigree data file (no need for the missing variables)*/
		if (fam != fam0){ /*They are not equal when have new family ID*/
			family[fam0].n_founder = n_founder;
			family[fam0].n_non_founder = n_non_founder;
			n_person = n_founder + n_non_founder;
			family[fam0].n_person = n_person;
			family[fam0].n_total = n_person;
			n_founder = n_non_founder = 0;

			/*Allocation of matrices and vectors for family ID fam0*/
			family[fam0].parent = imatrix(1, n_person, 0, 1);
			family[fam0].cov = dmatrix(1, n_person, 1, n_cov);
			family[fam0].genos = imatrix(1, n_person, 0, 1);
			family[fam0].phi = dmatrix(1, n_person, 1, n_person);
			family[fam0].delta_7 = dmatrix(1, n_person, 1, n_person);
			family[fam0].trait = dvector(1, n_person);
			family[fam0].mu_hat = dvector(1, n_person);
			family[fam0].y_perm = dvector(1, n_person);
					}
		fam0 = fam;

		if (fa == 0 && mo == 0) {
			n_founder++;
		}
		else if (fa != 0 && mo != 0) {
			n_non_founder++;
		}
		else {
			printf("ERROR: Father or mother of person %d is missing from pedigree file\n", ind);
		}
	}
	/*For the last family*/
	family[fam].n_founder = n_founder;
	family[fam].n_non_founder = n_non_founder;
	n_person = n_founder + n_non_founder;
	family[fam].n_person = n_person;
	family[fam].n_total = n_person;

	/*Allocation of matrices and vectors for family ID fam0*/
	family[fam].parent = imatrix(1, n_person, 0, 1);
	family[fam].cov = dmatrix(1, n_person, 1, n_cov);
	family[fam].genos = imatrix(1, n_person, 0, 1);
	family[fam].phi = dmatrix(1, n_person, 1, n_person);
	family[fam].delta_7 = dmatrix(1, n_person, 1, n_person);
	family[fam].trait = dvector(1, n_person);
	family[fam].mu_hat = dvector(1, n_person);
	family[fam].y_perm = dvector(1, n_person);

	fclose(pedfile);

} /* end pednum */

/* This function (readped) reads the pedigree information into **parent
   parent[person][which] - parent id of person, which = 0 for mother, which = 1 for father
   */

void readped(char* filename, struct FAMILY *family) {
	int i, j, k, fam, *ind, *fa, *mo, *sex, *age, n_person;

	n_person = family[1].n_person; /*Assumes all families have same number of individuals*/
	FILE *pedfile;

	ind = ivector(1, n_person);
	fa = ivector(1, n_person);
	mo = ivector(1, n_person);
	sex = ivector(1, n_person);
	age = ivector(1, n_person);

	pedfile = fopen(filename, "r");

	if (pedfile == NULL) {
		printf("Can't open pedigree file\n");
		exit(1);
	}


	k = 1;
	while (!feof(pedfile)) {
		fscanf(pedfile, "%d %d %d %d %d %d \n", &fam, &ind[k], &fa[k], &mo[k], &sex[k], &age[k]);

		if (k == n_person){
			k = 0;
			for (i = 1; i <= n_person; i++){
				family[fam].cov[i][1] = 1;
				family[fam].cov[i][2] = sex[i];
				family[fam].cov[i][3] = age[i];
			}

			ind[0] = 0;
			for (i = 1; i <= n_person; i++) {
				for (j = 0; j < i; j++) {
					if (fa[i] == ind[j]) { // Father
						family[fam].parent[i][0] = j;
					}
					if (mo[i] == ind[j]) { // Mother
						family[fam].parent[i][1] = j;
					}
				}
			}
		}
		k++;
	}

	free_ivector(ind, 1, n_person);
	free_ivector(fa, 1, n_person);
	free_ivector(mo, 1, n_person);
	free_ivector(sex, 1, n_person); /* Added by JM */
	free_ivector(age, 1, n_person); /* Added by JM */

	fclose(pedfile);
} /* end readped */

/* Read in the kinship and inbreeding coefficients from the output of KinInbCoef.C
   Store the kinship/inbreeding coefficient between pairs of individuals ind1 and ind2
   in a matrix storekin[ind1][ind2], this function assumes that everyone in the pedigree
   file is also in the kinship/inbreeding coefficient file
   */
void readkininb(char *filename, struct FAMILY *family) {

	FILE *kinfile;
	kinfile = fopen(filename, "r");

	if (kinfile == NULL) {
		printf("Can't open kinship and inbreeding coefficient file\n");
		exit(1);
	}

	int n_person;
	n_person = family[1].n_person;
	double kcoef = 0.0, dcoef = 0.0;
	int fam = 0, ind1 = 0, ind2 = 0, sub1 = 0, sub2 = 0;


	while (!feof(kinfile)) {
		if (fscanf(kinfile, "%d %d %d %lf %lf ", &fam, &ind1, &ind2, &kcoef, &dcoef) != 5) {
			break;
		}

		sub1 = ind1 - (fam - 1)*n_person;
		sub2 = ind2 - (fam - 1)*n_person;
		family[fam].phi[sub1][sub2] = kcoef;
		family[fam].phi[sub2][sub1] = kcoef;
		family[fam].delta_7[sub1][sub2] = dcoef;
		family[fam].delta_7[sub2][sub1] = dcoef;
	}

	fclose(kinfile);
} /* end readkininb */

// Store haplotype information from COSI to be use for founders' haplotypes in pedigree
void readhaplo(char *filename, int **haplotypes, int n_sites){
	int  i, j;
	int lines_in_haplo = 10000;

	FILE *haplofile;
	haplofile = fopen(filename, "r");
	if ((haplofile == NULL)) {
		printf("Can't open haplotype file\n");
		exit(1);
	}

	for (i = 1; i <= lines_in_haplo; i++){
		for (j = 1; j <= n_sites; j++){
			if (j == n_sites)
				fscanf(haplofile, "%d\n", &haplotypes[i][j]);
			else
				fscanf(haplofile, "%d ", &haplotypes[i][j]);
		}
	}
	fclose(haplofile);

}

