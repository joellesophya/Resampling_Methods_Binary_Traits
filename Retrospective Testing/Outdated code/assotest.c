#include "assotest.h"
#include "rand1.h"

// perm: indicates whether the response is the trait realization (perm=0) or permutation replicate (perm=1)
void VT(struct FAMILY *family, struct RV_INFO rv, int n_set, struct DATA_STRUCT data_struct, struct STATS stats, int n_perm, long *seed){

	int i, j, k, fam, m, uni, T, founder_allele, i_set, n0_perm = 1000;
	double tol = data_struct.tol, tmp = 0, cumsum_denum = 0, p_vals, tmp_p = 0;

	int n_person = family[1].n_person;
	int n_f = family[1].n_founder;

	int **genos = family[1].genos;
	int **haplo = rv.haplo;
	int **line_in_haplo = rv.lines_haplo;
	int **marknum = rv.poly_rv;

	int n_fam = data_struct.n_fam;
	int n_marker = data_struct.n_marker;
	int n_total = data_struct.n_total;
	int *n_case = data_struct.n_cases;
	double *sig_hat = data_struct.sig_hat;
	
	double *MAF = stats.MAF;
	double *sorted_MAF = stats.sorted_MAF;
	int *unique_MAF = stats.unique_MAF;
	double *threshold = stats.threshold;
	double *zVT = stats.zVT;
	double *cumsum_num = stats.cumsum_num;
	
	double ***G = d3tensor(1, n_fam, 1, n_person, 1, n_marker);
	
	for( i_set = 1; i_set <= n_set; i_set++){
		for (p_vals = 0, i = 1; i <= (n_perm + 1); i++) zVT[i] = 0, cumsum_num[i] = 0;
		
		// Step 1a: Compute MAF for each RV in the set
		for (i = 1; i <= n_marker; i++){
			threshold[i] = 0, m = 0;

			for (fam = 1; fam <= n_fam; fam++){
				genos = family[fam].genos;

				for (j = 1; j <= n_person; j++){
					founder_allele = genos[j][0] + (n_f * 2)*(fam - 1);
					m += (haplo[line_in_haplo[i_set][founder_allele]][marknum[i_set][i]] == 1) ? 1 : 0;
					founder_allele = genos[j][1] + (n_f * 2)*(fam - 1);
					m += (haplo[line_in_haplo[i_set][founder_allele]][marknum[i_set][i]] == 1) ? 1 : 0;
				}
			}
			MAF[i] = (1.0 + m) / (2.0 * (n_total + 1.0));
			sorted_MAF[i] = MAF[i];
		}

		// Step 1b: Get A where A=[G1,...,Gm]
		for (fam = 1; fam <= n_fam; fam++){
			genos = family[fam].genos;

			for (j = 1; j <= n_person; j++){
				for (i = 1; i <= n_marker; i++){
					tmp = 0;
					founder_allele = genos[j][0] + (n_f * 2)*(fam - 1);
					tmp += (haplo[line_in_haplo[i_set][founder_allele]][marknum[i_set][i]] == 1) ? 1 : 0;

					founder_allele = genos[j][1] + (n_f * 2)*(fam - 1);
					tmp += (haplo[line_in_haplo[i_set][founder_allele]][marknum[i_set][i]] == 1) ? 1 : 0;

					G[fam][j][i] = tmp;
				}
			}
		}

		// Step 2a: Determine the thresholds to use (sorted unique MAFs)
		qsort(&sorted_MAF[1], n_marker, sizeof(sorted_MAF[1]), cmp);

		for (i = 1, uni = 0; i <= n_marker; i++){
			if (i == 1) unique_MAF[i] = 1, uni++;
			else if ((sorted_MAF[i] - sorted_MAF[i - 1]) < tol) unique_MAF[i] = 0;
			else unique_MAF[i] = 1, uni++;
		}

		for (i = 1, T = 1; i <= n_marker; i++){
			if (unique_MAF[i] == 1)	threshold[T++] = sorted_MAF[i];
		}
		
		// ********  ADAPTIVE PERM PROCEDURE *************//
		// Step 2b: For each threshold value T, compute z(T):
		cumsum_denum = 0;
		for (T = 2; T <= uni; T++){
			for (i = 1; i <= n_marker; i++){
				if ((MAF[i] >= threshold[T - 1]) && (MAF[i] < threshold[T])){
					// marker MAF is not less than any of the previous thresholds
					//            but is less than the current threshold
					for(k = 1; k <= (n0_perm + 1); k++){
						// For each Y replicate, sum over sample (denominator is same for all Y)
						for (fam = 1; fam <= n_fam; fam++){
							for (j = 1; j <= n_person; j++){
								if (k == 1) cumsum_denum += G[fam][j][i] * G[fam][j][i];
								cumsum_num[k] += G[fam][j][i] * (family[fam].trait[k][j] - 1.0 * n_case[k] / n_total);
							}
						}
					}
				}
			}
			for(k = 1; k <= (n0_perm + 1); k++){
				tmp = fabs(cumsum_num[k]) / sqrt(cumsum_denum * sig_hat[k]);
				if (zVT[k] < tmp) zVT[k] = tmp;
			}
		}
		
		// Compute p-value from first n0_perm reps by comparing zVT between Y and Y replicates
		zVT[1] /= sqrt(n_case[1] * (n_total - n_case[1]) );
		for(k = 2; k <= (n0_perm + 1); k++){
			zVT[k] /= sqrt(n_case[k] * (n_total - n_case[k]) );
			if (zVT[1] < zVT[k]) p_vals++;
			else if (zVT[1] == zVT[k]) p_vals+= (rand1(seed) < .5) ? 1 : 0;
		}
		tmp_p = (p_vals + 1.0) / (n0_perm + 1.0);
		stats.pvals[i_set] = tmp_p;
		
		// Cutoff criterion based on Wald CI
		if ( (tmp_p - 1.96 * sqrt( tmp_p * (1 - tmp_p) / n0_perm ) ) < .1 ){
			// For p-value small enough, increase precision by adding more Y reps
			cumsum_denum = 0;
			for (T = 2; T <= uni; T++){
				for (i = 1; i <= n_marker; i++){
					if ((MAF[i] >= threshold[T - 1]) && (MAF[i] < threshold[T])){
						// marker MAF is not less than any of the previous thresholds
						//            but is less than the current threshold
						for( k = (n0_perm + 2); k <= (n_perm + 1) ; k++){
							// For each Y replicate, sum over sample (denominator is same for all Y)
							for (fam = 1; fam <= n_fam; fam++){
								for (j = 1; j <= n_person; j++){
									if (k == (n0_perm + 2)) cumsum_denum += G[fam][j][i] * G[fam][j][i];
									cumsum_num[k] += G[fam][j][i] * (family[fam].trait[k][j] - 1.0 * n_case[k] / n_total);
								}
							}
						}
					}
				}
				for( k = (n0_perm + 2); k <= (n_perm + 1); k++){
					tmp = fabs(cumsum_num[k]) / sqrt(cumsum_denum * sig_hat[k]);
					if (zVT[k] < tmp) zVT[k] = tmp;
				}
			}
			// Compute the p-value from additional reps by comparing zVT between Y and Y replicates
			for( k = (n0_perm + 2); k <= (n_perm + 1); k++){
				zVT[k] /= sqrt(n_case[k] * (n_total - n_case[k]) );
				if (zVT[1] < zVT[k]) p_vals++;
				else if (zVT[1] == zVT[k]) p_vals+= (rand1(seed) < .5) ? 1 : 0;
			}
			stats.pvals[i_set] = (p_vals + 1.0) / (n_perm + 1.0);
		}
	}
	free_d3tensor(G, 1, n_fam, 1, n_person, 1, n_marker + 1);	
}

// void wst(struct FAMILY *family, int **haplotypes, struct DATA_STRUCT data_struct, struct STATS *stats){
// int i, j, k, fam, m_u, n_u;
// double sum, x = 0;
// double tol = data_struct.tol;
// int n_person = family[1].n_person;
// int n_f = family[1].n_founder;
// int n_fam = data_struct.n_fam;
// int n_marker = data_struct.n_marker;
// int n_total = data_struct.n_total;

// double *trait = family[1].trait;
// int **genos = family[1].genos;

// // Step 1: For each variant, compute weight w_i
// double *w = dvector(1, n_marker);
// for (i = 1; i <= n_marker; i++){
// m_u = 0, n_u = 0;
// for (fam = 1; fam <= n_fam; fam++){
// genos = family[fam].genos;
// trait = family[fam].trait;
// for (j = 1; j <= n_person; j++){
// if (trait[j] < tol){
// m_u += (haplotypes[genos[j][0] + (n_f * 2)*(fam - 1)][i] == 1) ? 1 : 0;
// m_u += (haplotypes[genos[j][1] + (n_f * 2)*(fam - 1)][i] == 1) ? 1 : 0;
// n_u += 1;
// }
// }
// }
// w[i] = (1.0 + m_u) / (2.0 + 2.0 * n_u);
// w[i] = sqrt(n_total*w[i] * (1.0 - w[i]));
// }

// // Step 2: For each individual, compute score
// double **score = dmatrix(1, n_fam, 1, n_person);
// double *sorted_score = dvector(0, n_total - 1);
// k = 0;
// for (fam = 1; fam <= n_fam; fam++){
// genos = family[fam].genos;
// for (j = 1; j <= n_person; j++){
// score[fam][j] = 0.0;
// for (i = 1; i <= n_marker; i++){
// sum = 0.0;
// sum += (haplotypes[genos[j][0] + (n_f * 2)*(fam - 1)][i] == 1) ? 1.0 : 0.0;
// sum += (haplotypes[genos[j][1] + (n_f * 2)*(fam - 1)][i] == 1) ? 1.0 : 0.0;
// sum /= w[i];
// score[fam][j] += sum;
// }
// sorted_score[k] = score[fam][j];
// k++;
// }
// }

// // Step 3: Get the rank for all individuals and compute the sum of the ranks for cases
// qsort(sorted_score, n_total, sizeof(sorted_score[1]), cmp);
// double *ranks = dvector(1, n_total);

// i = 1;
// while (i <= n_total){
// k = 0;
// j = i;
// while (((sorted_score[j] - sorted_score[j - 1]) < tol) && (j<n_total)){
// k++, j++;
// }
// if (k > 0){
// n_u = 0;
// for (j = i; j <= (i + k); j++) n_u += j;
// sum = n_u / (k + 1.0);
// for (j = i; j <= (i + k); j++) ranks[j] = sum;
// i += k + 1;
// }
// else ranks[i] = i, i++;
// }

// for (fam = 1; fam <= n_fam; fam++){
// trait = family[fam].trait;
// for (j = 1; j <= n_person; j++){
// if (trait[j] > tol){
// for (k = 1; k <= n_total; k++){
// if (sorted_score[k - 1] == score[fam][j]){
// x += ranks[k];
// break;
// }
// }
// }
// }
// }

// stats->wss = x;

// free_dvector(w, 1, n_marker);
// free_dmatrix(score, 1, n_fam, 1, n_person);
// free_dvector(sorted_score, 0, n_total - 1);
// free_dvector(ranks, 1, n_total);
// }

int cmp(const void *a, const void *b){
	double aa = *(double*)a, bb = *(double*)b;
	return (aa > bb) - (aa < bb);
}
