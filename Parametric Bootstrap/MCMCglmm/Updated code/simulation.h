#ifndef SIMULATION_H
#define SIMULATION_H

#define MAXLEN 1024
#define MISSVAL -9.0
#define NUMTOL 1e-6
#define ANALYSIS_STATUS 1


struct FAMILY {
	int n_total;				/* total number individuals in family (SAME AS n_person) */
	int n_person;				/* number of phenotyped individuals in family */
	int n_founder, n_non_founder; /*Number of founders/non-founders*/

	int **parent;				//Matrix that contains info about number of parents: parent[i][j] = parent j ID for indiv i (j=0 for dad and 1 for mom)
	double **phi;				// phi[ind i][ind j] is 2*kinship/inbreeding matrix (i.e. phi[i][i]=1+h_i and phi[i][j]=2*kinship coef b/w i and j)
	double **delta_7;			/*Matrix that contains info about coefficient of identity*/
	double **cov;				/* cov[index among phenotyped][x] = cov value for ind for cov x */
									//column 2 - sex vector : sex[i] = gender of indiv i(= 1 for male and 2 for female)
									//column 3 - age vector : age[i] = age of indiv i
									//column 4 - normal covariate vector
	int **genos;                // genos[j][k]: founder haplotype # for individual j at allele k (k= 0/1)  so (n_person*2) matrix

	double **trait;				/* trait[Y-replicate][index among phenotyped] */
	double *y;					/* y[index] = 0.5x(no. of alleles 1);
								index = 1...n_geno_pheno, 1+n_geno_pheno...n_geno_notpheno+n_geno_pheno */
	double *mu_hat;             // Estimate of the mean of Y in GLMMM model

	////////////////////////////////////////////////////
	// For the permutation replicates
	double *y_perm; // Vector that holds Y permutation
};

struct SIM_PARAM {
	int n_sites;
	int n_marker;
	double add_var;
	double err_var;
	double dom_var;
	double sex_effect;
	double age_effect;
	int cov_type;
	int err_type;
	int fam0;
};

struct DATA_STRUCT {
	int n_fam;		/* total number of families */
	int n_marker;		/* total number of markers */
	int n_total;		/* total number of individuals (across all families) */
	int n_cov;		/* total number of covariates including intercept */
	int *n_cases;
	double tol;
	double **res_Y;
};

struct RV_INFO{
	int **haplo; // Stores haplotypes from COSI (10000*n_sites matrix)
	int **lines_haplo; // Stores row info for haplotypes for each set of RVs and each founder allele (n_set*2*n_founder_total matrix)
	int **poly_rv; // Stores info about polymorphic RV's for each set
	int n_set; // Number of RV sets being simulated (each set will correspond to a statistic)
};

struct STATS {
	double wss;
	double vt;
	double *mles;
	double *vt_obs;
	double *MAF;
	double *sorted_MAF;
	int *unique_MAF;
	double *threshold;
	double *zVT;
	double *cumsum_num;
	double *pvals;
};

struct PHI_SVD{ // Contains SVD decomp of family relatedness matrix (2*kinship coeff)
	double **u;
	double *lambda;
	double **sigma; // Estimate of omega = xi*Phi + (1-xi)*I
};

struct TMP_PERM{
	double *delta;
	double *delta_perm;
	double **v1;
	int nonzero;
	double **cholesky;
	double *v1Delta; // For obtaining the permuted phenotype vector
	double **v1Delta_f;
};

#endif
