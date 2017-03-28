#include "datacheck.h"
#include <string.h>
#include <stdlib.h>

void format_checks(char *pedfilename, struct DATA_STRUCT *data_struct) {

	/* define, initialize, and allocate memory for variables */
	int n_cov, i_cov, total, n_fam, fam_id0, fam_id, ind_id, fa, mo;
	int i, k;
	int ret;
	int sex;
	double cov_tmp = -9.0;
	double tol = data_struct->tol;
	int n_read = 0; /* to keep track of how much of the line has been read with sscanf */
	int n_read_tmp = 0; /* to store how much was read in the last call to sscanf */

	int col = 0;
	char car;

	/* check if files are present and can be opened input files */
	FILE *pedfile;
	pedfile = fopen(pedfilename, "r");
	if (pedfile == NULL) {
		printf("ERROR: Can't open pedigree and trait file\n");
		exit(1);
	}
	fclose(pedfile);

	/* read through the pedfile to count the number of columns, to then find the number of covariates */
	pedfile = fopen(pedfilename, "r");

	car = 0;
	while (car != EOF && car != '\n') {
		car = getc(pedfile);
		if (car == ' ' || car == '\t') {
			col++;
		}
	}
	col++;
	n_cov = col - 4 + 1; /* no. of covariates + 1 (+1 for the intercept which is assumed to be absent in pedfile) */
	fclose(pedfile);

	/* read through the pedfile to create lists of pedigree and individual ids,
	 count number of families (n_fam) and number of individuals (total) */
	fam_id0 = 0;// strdup("0");
	total = 0;
	n_fam = 0;
	pedfile = fopen(pedfilename, "r");
	while (!feof(pedfile)) {

		if (n_cov > 1){
			fscanf(pedfile, "%d %d %d %d %d ", &fam_id, &ind_id, &fa, &mo, &sex);
		}
		if (n_cov == 1) {
			fscanf(pedfile, "%d %d %d %d %d \n", &fam_id, &ind_id, &fa, &mo, &sex);
		}

		if (n_cov > 2) {
			for (i_cov = 3; i_cov < n_cov; i_cov++) {
				fscanf(pedfile, "%lf ", &cov_tmp);
			}
			fscanf(pedfile, "%lf \n", &cov_tmp);
		}
		total++;

		if (fam_id != fam_id0){
			n_fam++;
			fam_id0 = fam_id;
		}

	}

	fclose(pedfile);

	///* fill data_struct */
	data_struct->n_fam = n_fam;
	data_struct->n_total = total;
	data_struct->n_cov = n_cov;
}



