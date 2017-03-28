#include "assotest.h"

// perm: indicates whether the response is the trait realization (perm=0) or permutation replicate (perm=1)
// void VT(struct FAMILY *family, struct RV_INFO rv, int i_set, struct DATA_STRUCT data_struct, struct PHI_SVD phi_svd, struct STATS stats, int perm){

	// int i, j, k, fam, m, uni, T, g_ij, tmp_d, founder_allele;
	// double tol = data_struct.tol, tmp = 0, mean_trait = 0, cumsum_num = 0, cumsum_denum = 0, tmp_num = 0, tmp_denum = 0, sum_trait = 0.0;

	// double *trait;
	// int n_person = family[1].n_person;
	// int n_f = family[1].n_founder;

	// double **Sigma = phi_svd.sigma;
	// int **genos = family[1].genos;
	// int **haplo = rv.haplo;
	// int **line_in_haplo = rv.lines_haplo;
	// int **marknum = rv.poly_rv;

	// int n_fam = data_struct.n_fam;
	// int n_marker = data_struct.n_marker;
	// int n_total = data_struct.n_total;
	
	// double *MAF = stats.MAF;
	// double *sorted_MAF = stats.sorted_MAF;
	// int *unique_MAF = stats.unique_MAF;
	// double *threshold = stats.threshold;
	// double *z = stats.z;
	// double *vt = stats.obs_stat;

	// // Step 1a: Compute MAF for each variant
	// sum_trait = 0.0;
	// for (i = 1; i <= n_marker; i++){
		// threshold[i] = 0;
		// z[i] = 0;
		// m = 0;

		// for (fam = 1; fam <= n_fam; fam++){
			// genos = family[fam].genos;

			// if (perm == 0) trait = family[fam].trait;
			// else  trait = family[fam].y_perm;

			// for (j = 1; j <= n_person; j++){
				// founder_allele = genos[j][0] + (n_f * 2)*(fam - 1);
				// m += (haplo[line_in_haplo[i_set][founder_allele]][marknum[i_set][i]] == 1) ? 1 : 0;
				// founder_allele = genos[j][1] + (n_f * 2)*(fam - 1);
				// m += (haplo[line_in_haplo[i_set][founder_allele]][marknum[i_set][i]] == 1) ? 1 : 0;
				// if (i == 1) sum_trait += trait[j]; // To compute Y_bar
			// }
		// }

		// MAF[i] = (1.0 + m) / (2.0 * (n_total + 1.0));
		// sorted_MAF[i] = MAF[i];
	// }
	// mean_trait = sum_trait / n_total*1.0; // Mean trait value (if binary trait then = proportion of cases)

	// // Step 1b: Get C^t*A where A=[G1,...,Gm,Y-Ybar] for each family and Phi
	// double ***GY = d3tensor(1, n_fam, 1, n_person, 1, n_marker + 1);
	// double **chol = dmatrix(1, n_person, 1, n_person);
	// for (fam = 1; fam <= n_fam; fam++){
		// genos = family[fam].genos;
		// if (perm == 0) trait = family[fam].trait;
		// else  trait = family[fam].y_perm;

		// for (j = 1; j <= n_person; j++){
			// for (i = 1; i <= n_marker; i++){
				// tmp = 0;
				// founder_allele = genos[j][0] + (n_f * 2)*(fam - 1);
				// tmp += (haplo[line_in_haplo[i_set][founder_allele]][marknum[i_set][i]] == 1) ? 1 : 0;

				// founder_allele = genos[j][1] + (n_f * 2)*(fam - 1);
				// tmp += (haplo[line_in_haplo[i_set][founder_allele]][marknum[i_set][i]] == 1) ? 1 : 0;

				// GY[fam][j][i] = tmp;
			// }
			// GY[fam][j][n_marker + 1] = trait[j] - mean_trait;
		// }
		// cholesky(Sigma, n_person, GY[fam], n_marker + 1, chol, GY[fam], 1);
	// }


	// // Step 2a: Determine the thresholds to use (sorted unique MAFs)
	// qsort(&sorted_MAF[1], n_marker, sizeof(sorted_MAF[1]), cmp);

	// for (i = 1, uni = 0; i <= n_marker; i++){
		// if (i == 1) unique_MAF[i] = 1, uni++;
		// else if ((sorted_MAF[i] - sorted_MAF[i - 1]) < tol) unique_MAF[i] = 0;
		// else unique_MAF[i] = 1, uni++;
	// }

	// for (i = 1, T = 1; i <= n_marker; i++){
		// if (unique_MAF[i] == 1)	threshold[T++] = sorted_MAF[i];
	// }

	// // Step 2b: For each threshold value T, compute z(T):
	// cumsum_num = cumsum_denum = 0.0;
	// // No need to compute z(T) for lowest threshold
	// for (T = 2; T <= uni; T++){
		// for (i = 1; i <= n_marker; i++){
			// if ((MAF[i] >= threshold[T - 1]) && (MAF[i] < threshold[T])){
				// // marker MAF is not less than any of the previous thresholds
				// //            but is less than the current threshold
				// tmp_num = tmp_denum = 0.0;
				// for (fam = 1; fam <= n_fam; fam++){
					// for (j = 1; j <= n_person; j++){
						// tmp_denum += GY[fam][j][i] * GY[fam][j][i];
						// tmp = GY[fam][j][i] * GY[fam][j][n_marker + 1];
						// tmp_num += tmp;
					// }
				// }
				// cumsum_num += tmp_num;
				// cumsum_denum += tmp_denum;
			// }
		// }
		// z[T - 1] = fabs(cumsum_num) / sqrt(cumsum_denum);
	// }

	// // Step 3: Compute z_max	
	// qsort(&z[1], uni - 1, sizeof(z[1]), cmp);
	
	// if (perm == 0) vt[i_set] = z[uni - 1];
	// z[1] = z[uni - 1];

	
	// free_d3tensor(GY, 1, n_fam, 1, n_person, 1, n_marker + 1);
	// free_dmatrix(chol, 1, n_person, 1, n_person);
	// //FILE *outfile;
	// //outfile = fopen("test.txt", "w");
	// //// Write file that is [Y G] (trait and genotype matrix)
	// //for (fam = 1; fam <= n_fam; fam++){
	// //	trait = family[fam].y_perm;
	// //	genos = family[fam].genos;
	// //	for (j = 1; j <= n_person; j++){
	// //		fprintf(outfile, "%lf\t", trait[j]);
	// //		for (i = 1; i <= n_marker; i++){
	// //			tmp_d = 0;
	// //			founder_allele = genos[j][0] + (n_f * 2)*(fam - 1);
	// //			tmp_d += (haplo[line_in_haplo[i_set][founder_allele]][marknum[i_set][i]] == 1) ? 1 : 0;
	// //			founder_allele = genos[j][1] + (n_f * 2)*(fam - 1);
	// //			tmp_d += (haplo[line_in_haplo[i_set][founder_allele]][marknum[i_set][i]] == 1) ? 1 : 0;
	// //			fprintf(outfile, "%d\t", tmp_d);
	// //		}
	// //		fprintf(outfile, "\n");
	// //	}
	// //}
	//	// Prints out [Y X] (trait and covariates (w/o intercept))
	//for (fam = 1; fam <= n_fam; fam++){
	//	trait = family[fam].trait;
	//	for (j = 1; j <= n_person; j++){
	//		fprintf(outfile, "%lf\t", trait[j]);
	//		fprintf(outfile, "%d\t%d\n", family[fam].age[j], family[fam].sex[j]);
	//	}
	//}
	//fclose(outfile);
// }

void wst(struct FAMILY *family, struct RV_INFO rv, int i_set, struct DATA_STRUCT data_struct, struct STATS stats, int perm){
	int i, j, k, fam, m_u, n_u, founder_allele, n_match;
	double sum, x, sum_ranks;
	double tol = data_struct.tol;
	
	int n_person = family[1].n_person;
	int n_f = family[1].n_founder;
	int n_fam = data_struct.n_fam;
	int n_marker = data_struct.n_marker;
	int n_total = data_struct.n_total;

	double *trait = family[1].trait;
	int **genos = family[1].genos;
	
	int **haplo = rv.haplo;
	int **line_in_haplo = rv.lines_haplo;
	int **marknum = rv.poly_rv;
	
	double *w = stats.z;	
	double *wst = stats.obs_stat;
	double **score = stats.score;
	double *sorted_score = stats.sorted_score;
	double *score_ranks = stats.score_ranks;
	
	// Step 1: For each variant, compute weight w_i
	for (i = 1; i <= n_marker; i++){
		m_u = 0, n_u = 0;
		for (fam = 1; fam <= n_fam; fam++){
			genos = family[fam].genos;
			trait = (perm == 0)? family[fam].trait : family[fam].y_perm;
			
			for (j = 1; j <= n_person; j++){
				if (trait[j] < tol){			
					founder_allele = genos[j][0] + (n_f * 2)*(fam - 1);
					m_u += (haplo[line_in_haplo[i_set][founder_allele]][marknum[i_set][i]] == 1) ? 1 : 0;
					founder_allele = genos[j][1] + (n_f * 2)*(fam - 1);
					m_u += (haplo[line_in_haplo[i_set][founder_allele]][marknum[i_set][i]] == 1) ? 1 : 0;
					n_u ++;
				}
			}
		}
		w[i] = (1.0 + m_u) / (2.0 * (1 + n_u));
		w[i] = sqrt(n_total * w[i] * (1.0 - w[i]));
	}
	
	// Step 2: For each individual, compute score	
	for (fam = 1, k = 1; fam <= n_fam; fam++){
		genos = family[fam].genos;
		for (j = 1; j <= n_person; j++, k++){
			for (i = 1, score[fam][j] = 0; i <= n_marker; i++){	
				sum = 0.0;			
				founder_allele = genos[j][0] + (n_f * 2)*(fam - 1);
				sum += (haplo[line_in_haplo[i_set][founder_allele]][marknum[i_set][i]] == 1) ? 1 : 0;
				founder_allele = genos[j][1] + (n_f * 2)*(fam - 1);
				sum += (haplo[line_in_haplo[i_set][founder_allele]][marknum[i_set][i]] == 1) ? 1 : 0;
				sum /= w[i];
				score[fam][j] += sum;
			}
			sorted_score[k] = score[fam][j];
			score_ranks[k] = k * 1.0;
		}
	}

	// Step 3: Get the rank for all individuals and compute the sum of the ranks for cases
	qsort(&sorted_score[1], n_total, sizeof(sorted_score[1]), cmp);
	
	// Modify ranks to account for matches by taking average of ranks
	for(i = 1; i <= n_total ; i++){
		n_match = 0, sum_ranks = score_ranks[i], k = i;
		if(i < n_total){
			while(fabs(sorted_score[i] - sorted_score[i+1]) < tol){
				n_match++;
				sum_ranks += score_ranks[++i];
				if (i == n_total) break;
			}
		}
		for(j = 0; j <= n_match ; j++) score_ranks[k + j] = sum_ranks / (n_match + 1.0);
	}

	for (fam = 1, x = 0; fam <= n_fam; fam++){
		trait = (perm == 0)? family[fam].trait : family[fam].y_perm;
		
		for (j = 1; j <= n_person; j++){
			if (trait[j] > tol){ // For cases
				for (k = 1; k <= n_total; k++){
					if (sorted_score[k] == score[fam][j]){
						x += score_ranks[k];
						break;  // Only need to find first match
					}
				}
			}
		}
	}
	
	if (perm == 0) wst[i_set] = x;
	w[1] = x; // This is the WST statistic
	
	// outfile = fopen("test.txt", "w");
	// // Write file that is [Y G] (trait and genotype matrix)
	// for (fam = 1; fam <= n_fam; fam++){
		// trait = (perm == 0)? family[fam].trait : family[fam].y_perm;
		// genos = family[fam].genos;
		// for (j = 1; j <= n_person; j++){
			// fprintf(outfile, "%f\t", trait[j]);
			// for (i = 1; i <= n_marker; i++){
				// m_u = 0;
				// founder_allele = genos[j][0] + (n_f * 2)*(fam - 1);
				// m_u += (haplo[line_in_haplo[i_set][founder_allele]][marknum[i_set][i]] == 1) ? 1 : 0;
				// founder_allele = genos[j][1] + (n_f * 2)*(fam - 1);
				// m_u += (haplo[line_in_haplo[i_set][founder_allele]][marknum[i_set][i]] == 1) ? 1 : 0;
				// fprintf(outfile, "%d\t", m_u);
			// }
			// fprintf(outfile, "\n");
		// }
	// }
		// // Prints out [Y X] (trait and covariates (w/o intercept))
	// for (fam = 1; fam <= n_fam; fam++){
		// trait = family[fam].trait;
		// for (j = 1; j <= n_person; j++){
			// fprintf(outfile, "%lf\t", trait[j]);
			// fprintf(outfile, "%d\t%d\n", family[fam].age[j], family[fam].sex[j]);
		// }
	// }
	// fclose(outfile);
	// exit(1);
}

int cmp(const void *a, const void *b){
	double aa = *(double*)a, bb = *(double*)b;
	return (aa > bb) - (aa < bb);
}
