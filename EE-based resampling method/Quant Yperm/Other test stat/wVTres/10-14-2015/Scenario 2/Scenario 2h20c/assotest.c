#include "assotest.h"

// perm: indicates whether the response is the trait realization (perm=0) or permutation replicate (perm=1)
void VT(struct FAMILY *family, struct RV_INFO rv, int i_set, struct DATA_STRUCT data_struct, struct PHI_SVD phi_svd, struct STATS *stats,
	double *MAF, double *sorted_MAF, int *unique_MAF, double *threshold, double *z, int perm){

	int i, j, k, fam, m, uni, T, g_ij, tmp_d, founder_allele;
	double tol = data_struct.tol, tmp = 0, cumsum_num = 0, cumsum_denum = 0, tmp_num = 0, tmp_denum = 0;

	double *trait;
	int n_person = family[1].n_person;
	int n_f = family[1].n_founder;

	int **genos = family[1].genos;
	int **haplo = rv.haplo;
	int **line_in_haplo = rv.lines_haplo;
	int **marknum = rv.poly_rv;

	int n_fam = data_struct.n_fam;
	int n_cov = data_struct.n_cov;
	int n_marker = data_struct.n_marker;
	int n_total = data_struct.n_total;
	double **resid= data_struct.res_Y;
	double *bhat = data_struct.beta_hat;

	// Step 1a: Compute MAF for each variant
	threshold[0] = 0;
	for (i = 1; i <= n_marker; i++){
		threshold[i] = 0;
		z[i - 1] = 0;
		m = 0;

		for (fam = 1; fam <= n_fam; fam++){
			genos = family[fam].genos;

			if (perm == 0) trait = family[fam].trait;
			else  trait = family[fam].y_perm;

			for (j = 1; j <= n_person; j++){
				founder_allele = genos[j][0] + (n_f * 2)*(fam - 1);
				m += (haplo[line_in_haplo[i_set][founder_allele]][marknum[i_set][i]] == 1) ? 1 : 0;
				founder_allele = genos[j][1] + (n_f * 2)*(fam - 1);
				m += (haplo[line_in_haplo[i_set][founder_allele]][marknum[i_set][i]] == 1) ? 1 : 0;
			}
		}
		MAF[i] = (1.0 + m) / (2.0 * (n_total + 1.0));
		sorted_MAF[i - 1] = MAF[i];
	}

	// Step 1b: Get A where A=[G1,...,Gm,Y-Yhat] for each family 
	double ***GY = d3tensor(1, n_fam, 1, n_person, 1, n_marker + 1);
	for (fam = 1; fam <= n_fam; fam++){
		genos = family[fam].genos;
		if (perm == 0) trait = family[fam].trait;
		else  trait = family[fam].y_perm;

		for (j = 1; j <= n_person; j++){
			for (i = 1; i <= n_marker; i++){
				tmp = 0;
				founder_allele = genos[j][0] + (n_f * 2)*(fam - 1);
				tmp += (haplo[line_in_haplo[i_set][founder_allele]][marknum[i_set][i]] == 1) ? 1 : 0;

				founder_allele = genos[j][1] + (n_f * 2)*(fam - 1);
				tmp += (haplo[line_in_haplo[i_set][founder_allele]][marknum[i_set][i]] == 1) ? 1 : 0;

				GY[fam][j][i] = tmp;
			}
			GY[fam][j][n_marker + 1] = resid[fam][j];
		}
	}


	// Step 2a: Determine the thresholds to use (sorted unique MAFs)
	qsort(sorted_MAF, n_marker, sizeof(sorted_MAF[1]), cmp);

	uni = 0;
	for (i = 1; i <= n_marker; i++){
		if (i == 1) unique_MAF[i] = 1, uni++;
		else if ((sorted_MAF[i - 1] - sorted_MAF[i - 2]) < tol) unique_MAF[i] = 0;
		else unique_MAF[i] = 1, uni++;
	}

	//double *threshold = dvector(1, uni);
	T = 1;
	for (i = 1; i <= n_marker; i++){
		if (unique_MAF[i] == 1){
			threshold[T] = sorted_MAF[i - 1];
			T++;
		}
	}

	// Step 2b: For each threshold value T, compute z(T):
	cumsum_num = cumsum_denum = 0.0;
	// No need to compute z(T) for lowest threshold
	for (T = 2; T <= uni; T++){
		for (i = 1; i <= n_marker; i++){
			if ((MAF[i] >= threshold[T - 1]) && (MAF[i] < threshold[T])){
				// marker MAF is not less than any of the previous thresholds
				//            but is less than the current threshold
				tmp_num = tmp_denum = 0.0;
				for (fam = 1; fam <= n_fam; fam++){
					for (j = 1; j <= n_person; j++){

						tmp_denum += GY[fam][j][i] * GY[fam][j][i];
						tmp = GY[fam][j][i] * GY[fam][j][n_marker + 1];
						tmp_num += tmp;
					}
				}
				cumsum_num += tmp_num;
				cumsum_denum += tmp_denum;
			}
		}
		z[T - 2] = fabs(cumsum_num) / sqrt(cumsum_denum * bhat[n_cov + 1]);
	}

	// Step 3: Compute z_max	
	qsort(z, uni - 1, sizeof(z[1]), cmp);
	stats->vt = z[uni - 2];

	free_d3tensor(GY, 1, n_fam, 1, n_person, 1, n_marker + 1);
	
	//FILE *outfile;
	//outfile = fopen("test.txt", "w");
	//// Write file that is [Y G] (trait and genotype matrix)
	//for (fam = 1; fam <= n_fam; fam++){
	//	trait = family[fam].y_perm;
	//	genos = family[fam].genos;
	//	for (j = 1; j <= n_person; j++){
	//		fprintf(outfile, "%lf\t", trait[j]);
	//		for (i = 1; i <= n_marker; i++){
	//			tmp_d = 0;
	//			founder_allele = genos[j][0] + (n_f * 2)*(fam - 1);
	//			tmp_d += (haplo[line_in_haplo[i_set][founder_allele]][marknum[i_set][i]] == 1) ? 1 : 0;
	//			founder_allele = genos[j][1] + (n_f * 2)*(fam - 1);
	//			tmp_d += (haplo[line_in_haplo[i_set][founder_allele]][marknum[i_set][i]] == 1) ? 1 : 0;
	//			fprintf(outfile, "%d\t", tmp_d);
	//		}
	//		fprintf(outfile, "\n");
	//	}
	//}
	//	// Prints out [Y X] (trait and covariates (w/o intercept))
	//for (fam = 1; fam <= n_fam; fam++){
	//	trait = family[fam].trait;
	//	for (j = 1; j <= n_person; j++){
	//		fprintf(outfile, "%lf\t", trait[j]);
	//		fprintf(outfile, "%d\t%d\n", family[fam].age[j], family[fam].sex[j]);
	//	}
	//}
	//fclose(outfile);
}

void wst(struct FAMILY *family, int **haplotypes, struct DATA_STRUCT data_struct, struct STATS *stats){
	int i, j, k, fam, m_u, n_u;
	double sum, x = 0;
	double tol = data_struct.tol;
	int n_person = family[1].n_person;
	int n_f = family[1].n_founder;
	int n_fam = data_struct.n_fam;
	int n_marker = data_struct.n_marker;
	int n_total = data_struct.n_total;

	double *trait = family[1].trait;
	int **genos = family[1].genos;

	// Step 1: For each variant, compute weight w_i
	double *w = dvector(1, n_marker);
	for (i = 1; i <= n_marker; i++){
		m_u = 0, n_u = 0;
		for (fam = 1; fam <= n_fam; fam++){
			genos = family[fam].genos;
			trait = family[fam].trait;
			for (j = 1; j <= n_person; j++){
				if (trait[j] < tol){
					m_u += (haplotypes[genos[j][0] + (n_f * 2)*(fam - 1)][i] == 1) ? 1 : 0;
					m_u += (haplotypes[genos[j][1] + (n_f * 2)*(fam - 1)][i] == 1) ? 1 : 0;
					n_u += 1;
				}
			}
		}
		w[i] = (1.0 + m_u) / (2.0 + 2.0 * n_u);
		w[i] = sqrt(n_total*w[i] * (1.0 - w[i]));
	}

	// Step 2: For each individual, compute score
	double **score = dmatrix(1, n_fam, 1, n_person);
	double *sorted_score = dvector(0, n_total - 1);
	k = 0;
	for (fam = 1; fam <= n_fam; fam++){
		genos = family[fam].genos;
		for (j = 1; j <= n_person; j++){
			score[fam][j] = 0.0;
			for (i = 1; i <= n_marker; i++){
				sum = 0.0;
				sum += (haplotypes[genos[j][0] + (n_f * 2)*(fam - 1)][i] == 1) ? 1.0 : 0.0;
				sum += (haplotypes[genos[j][1] + (n_f * 2)*(fam - 1)][i] == 1) ? 1.0 : 0.0;
				sum /= w[i];
				score[fam][j] += sum;
			}
			sorted_score[k] = score[fam][j];
			k++;
		}
	}

	// Step 3: Get the rank for all individuals and compute the sum of the ranks for cases
	qsort(sorted_score, n_total, sizeof(sorted_score[1]), cmp);
	double *ranks = dvector(1, n_total);

	i = 1;
	while (i <= n_total){
		k = 0;
		j = i;
		while (((sorted_score[j] - sorted_score[j - 1]) < tol) && (j<n_total)){
			k++, j++;
		}
		if (k > 0){
			n_u = 0;
			for (j = i; j <= (i + k); j++) n_u += j;
			sum = n_u / (k + 1.0);
			for (j = i; j <= (i + k); j++) ranks[j] = sum;
			i += k + 1;
		}
		else ranks[i] = i, i++;
	}

	for (fam = 1; fam <= n_fam; fam++){
		trait = family[fam].trait;
		for (j = 1; j <= n_person; j++){
			if (trait[j] > tol){
				for (k = 1; k <= n_total; k++){
					if (sorted_score[k - 1] == score[fam][j]){
						x += ranks[k];
						break;
					}
				}
			}
		}
	}

	stats->wss = x;

	free_dvector(w, 1, n_marker);
	free_dmatrix(score, 1, n_fam, 1, n_person);
	free_dvector(sorted_score, 0, n_total - 1);
	free_dvector(ranks, 1, n_total);
}

int cmp(const void *a, const void *b){
	double aa = *(double*)a, bb = *(double*)b;
	return (aa > bb) - (aa < bb);
}