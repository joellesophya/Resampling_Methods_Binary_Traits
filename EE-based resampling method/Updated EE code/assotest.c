#include "assotest.h"
#include "rand1.h"
double *array;

void VT(struct FAMILY *family, struct RV_INFO rv, int n_set, struct DATA_STRUCT data_struct, double *pvals, int n_perm, long *seed){

	int i, j, k, fam, i_set, m, uni, T, founder_allele, p_vals, n0_perm = n_perm / 10, **genos;
	double tol = data_struct.tol, tmp, cumsum_denum, cutoff = .1, Qz = 2.58, CI_lb, varY; //99% conf. lvl

	int n_person = family[1].n_person;
	int n_f = family[1].n_founder;

	int **haplo = rv.haplo;
	int **line_in_haplo = rv.lines_haplo;
	int **marknum = rv.poly_rv;

	int n_fam = data_struct.n_fam;
	int n_marker = data_struct.n_marker;
	int n_total = data_struct.n_total;
	double *sumsY = data_struct.sumsY;
	
	double *MAF = dvector(1, n_marker);
	double *sorted_MAF = dvector(1, n_marker + 1);
	int *unique_MAF = ivector(1, n_marker + 1);
	double *threshold = dvector(1, n_marker + 1); // T = 1 has Tmax
	double *zVT = dvector(1, n_perm + 1);
	double *cumsum_num = dvector(1, n_perm + 1);
	double ***G = d3tensor(1, n_fam, 1, n_person, 1, n_marker);
	
	for( i_set = 1; i_set <= n_set; i_set++){
		// Step 1a: Get A where A=[G1,...,Gm] and Gi is in {0,1,2}
		for (fam = 1; fam <= n_fam; fam++){
			genos = family[fam].genos;

			for (j = 1; j <= n_person; j++){
				for (i = 1; i <= n_marker; i++){
					tmp = 0;
					founder_allele = genos[j][0] + (n_f * 2) * (fam - 1);
					tmp += (haplo[line_in_haplo[i_set][founder_allele]][marknum[i_set][i]] == 1) ? 1 : 0;

					founder_allele = genos[j][1] + (n_f * 2)*(fam - 1);
					tmp += (haplo[line_in_haplo[i_set][founder_allele]][marknum[i_set][i]] == 1) ? 1 : 0;

					G[fam][j][i] = tmp;
				}
			}
		}
		
		// Step 1b: Compute MAF for each RV in the set
		for (i = 1; i <= n_marker; i++){
			for ( threshold[i] = 0, m = 0, fam = 1; fam <= n_fam; fam++){
				for (j = 1; j <= n_person; j++) m += G[fam][j][i];
			}
			MAF[i] = (1.0 + m) / (2.0 * (n_total + 1));
			sorted_MAF[i] = MAF[i];
		}
		sorted_MAF[n_marker + 1] = 1;

		// Step 2a: Determine the thresholds to use (sorted unique MAFs)
		qsort(&sorted_MAF[1], n_marker + 1, sizeof(sorted_MAF[1]), cmp);

		for (uni = 0, i = 1; i <= (n_marker + 1); i++){
			if (i == 1) unique_MAF[i] = 1, uni++;
			else if ((sorted_MAF[i] - sorted_MAF[i - 1]) < tol) unique_MAF[i] = 0;
			else unique_MAF[i] = 1, uni++;
		}

		for (T = 1, i = 1; i <= (n_marker + 1); i++){
			if (unique_MAF[i] == 1) threshold[T++] = sorted_MAF[i];
		}
		
		// ********  ADAPTIVE PERM PROCEDURE *************//
		for (i = 1; i <= (n_perm + 1); i++) zVT[i] = cumsum_num[i] = 0; // Set num. to 0 FOR ALL n_perm replicates

		// Step 2b: For each threshold value T, compute z(T):
		for (cumsum_denum = 0, T = 2; T <= uni; T++){
			for (i = 1; i <= n_marker; i++){
				if ((MAF[i] >= (threshold[T - 1] - tol)) && (MAF[i] < (threshold[T] - tol))){
					// marker MAF is g.e.q. than the previous threshold but is less than the current one
					for(k = 1; k <= (n0_perm + 1); k++){
						// For each Y replicate, take G^T(Y-Ybar) & denominator is same for all Y
						for (fam = 1; fam <= n_fam; fam++){
							for (j = 1; j <= n_person; j++){
								if (k == 1) cumsum_denum += G[fam][j][i] * G[fam][j][i];
								cumsum_num[k] += G[fam][j][i] * (family[fam].trait[k][j] -  (sumsY[k] / n_total));
							}
						}
					}
				}
			}
			for(k = 1; k <= (n0_perm + 1); k++){
				varY = (sumsY[k + (n_perm + 1)] - sumsY[k] * sumsY[k] / n_total) / (n_total - 1);
				tmp = fabs(cumsum_num[k]) / sqrt(cumsum_denum * varY);
				if ((zVT[k] + tol) < tmp) zVT[k] = tmp;
			}
		}
		
		// Compare VT stat between Y and first n0_perm Y replicates to get p-value
		for(p_vals = 0, k = 2; k <= (n0_perm + 1); k++){
			if ((zVT[1] + tol) < zVT[k] ) p_vals++;
			else if (fabs(zVT[1] - zVT[k]) <= tol) p_vals += (rand1(seed) < .5) ? 1 : 0;
		}
		pvals[i_set] = (p_vals + 1.0) / (n0_perm + 1.0);
		
		// Cutoff criterion based on Normal Approx CI (not valid for p close to 0/1)
		CI_lb = pvals[i_set] - Qz * sqrt(pvals[i_set] * (1 - pvals[i_set]) / n0_perm);
		if (CI_lb <= (cutoff + tol)){
			// For p-values small enough, increase precision by adding more Y replicates
			for (cumsum_denum = 0, T = 2; T <= uni; T++){
				for (i = 1; i <= n_marker; i++){
					if ((MAF[i] >= (threshold[T - 1] - tol)) && (MAF[i] < (threshold[T] - tol))){
						// marker MAF is g.e.q. than the previous threshold but is less than the current one
						for ( k = (n0_perm + 2); k <= (n_perm + 1) ; k++){
							// For each Y replicate, take G^T(Y-Ybar) & denominator is same for all Y
							for (fam = 1; fam <= n_fam; fam++){
								for (j = 1; j <= n_person; j++){
									if (k == (n0_perm + 2)) cumsum_denum += G[fam][j][i] * G[fam][j][i];
									cumsum_num[k] += G[fam][j][i] * (family[fam].trait[k][j] - (sumsY[k] / n_total));
								}
							}
						}
					}
				}
				for( k = (n0_perm + 2); k <= (n_perm + 1); k++){
					varY = (sumsY[k + (n_perm + 1)] - sumsY[k] * sumsY[k] / n_total) / (n_total - 1);
					tmp = fabs(cumsum_num[k]) / sqrt(cumsum_denum * varY);
					if ((zVT[k] + tol) < tmp) zVT[k] = tmp;
				}
			}
			// Compute the p-value from additional replicates by comparing zVT between Y and Y replicates
			for( k = (n0_perm + 2); k <= (n_perm + 1); k++){
				if ((zVT[1] + tol) < zVT[k] ) p_vals++;
				else if (fabs(zVT[1]-zVT[k]) <= tol) p_vals+= (rand1(seed) < .5) ? 1 : 0;
			}
			pvals[i_set] = (p_vals + 1.0) / (n_perm + 1.0);
		}
	}
	
	free_dvector(threshold, 1, n_marker + 1);
	free_dvector(MAF, 1, n_marker);
	free_dvector(sorted_MAF, 1, n_marker + 1);
	free_ivector(unique_MAF, 1, n_marker + 1);
	free_dvector(zVT, 1, n_perm + 1);
	free_dvector(cumsum_num, 1, n_perm + 1);
	free_d3tensor(G, 1, n_fam, 1, n_person, 1, n_marker + 1);	
}

void ws(struct FAMILY *family, struct RV_INFO rv, int n_set, struct DATA_STRUCT data_struct, double *pvals, int n_perm, long *seed){
	
	int i, j, k, fam, m_u, n_u, founder_allele, i_y, i_set, n0_perm = n_perm / 10, **genos;
	double Gi, zWS, x, tol = data_struct.tol, p_vals, CI_lb, cutoff = .1, Qz = 2.58;
	
	int n_person = family[1].n_person;
	int n_f = family[1].n_founder;
	
	int n_fam = data_struct.n_fam;
	int n_marker = data_struct.n_marker;
	int n_total = data_struct.n_total;

	int **haplo = rv.haplo;
	int **line_in_haplo = rv.lines_haplo;
	int **marknum = rv.poly_rv;
	
	double *w = dvector(1, n_marker);
	double *score = dvector(1, n_total);
	int *score_ind = ivector(1, n_total);
	double *ranks = dvector(1, n_total);

	
	// For each RV set
	for( i_set = 1; i_set <= n_set; i_set++){
		// For each first n0_perm replicates (excluding obs Y)
		for( p_vals = zWS = 0, i_y = 1; i_y <= (n0_perm + 1); i_y++){
			
			// Step 1: For each variant, compute weight w_i
			for (i = 1; i <= n_marker; i++){
				for (m_u = 0, n_u = 0, fam = 1; fam <= n_fam; fam++){
					genos = family[fam].genos;

					for (j = 1; j <= n_person; j++){
						if (family[fam].trait[i_y][j] < tol){
							founder_allele = genos[j][0] + (n_f * 2)*(fam - 1);
							m_u += (haplo[line_in_haplo[i_set][founder_allele]][marknum[i_set][i]] == 1) ? 1 : 0;
							founder_allele = genos[j][1] + (n_f * 2)*(fam - 1);
							m_u += (haplo[line_in_haplo[i_set][founder_allele]][marknum[i_set][i]] == 1) ? 1 : 0;
							n_u++;
						}
					}
				}
				w[i] = (1.0 + m_u) / (2.0 * (1 + n_u));
				w[i] = sqrt(n_total * w[i] * (1 - w[i]));
			}
			
			// Step 2: For each individual, compute score
			for (k = 1, fam = 1; fam <= n_fam; fam++){
				genos = family[fam].genos;

				for (j = 1; j <= n_person; j++, k++){
					for (x = 0, i = 1; i <= n_marker; i++){
						founder_allele = genos[j][0] + (n_f * 2)*(fam - 1);
						Gi = (haplo[line_in_haplo[i_set][founder_allele]][marknum[i_set][i]] == 1) ? 1 : 0;
						founder_allele = genos[j][1] + (n_f * 2)*(fam - 1);
						Gi += (haplo[line_in_haplo[i_set][founder_allele]][marknum[i_set][i]] == 1) ? 1 : 0;
						Gi /= w[i];
						x += Gi;
					}
					score[k] = x, score_ind[k] = k;
				}
			}

			// Step 3: Get the rank for all individuals and compute the sum of the ranks for cases
			getRank(score, score_ind, ranks, n_total);
			for (x = 0, k = 1, fam = 1; fam <= n_fam; fam++){
				for (j = 1; j <= n_person; j++, k++){
					x += ranks[k] * family[fam].trait[i_y][j];
				}
			}
			if(i_y == 1) zWS = x;
			else{
				// Compare WS between Y and Y rep
				if ((zWS + tol) < x) p_vals++;
				else if (fabs(zWS - x) <= tol) p_vals += (rand1(seed) < .5) ? 1 : 0;
			}
		}
		pvals[i_set] = (p_vals + 1.0) / (n0_perm + 1.0);
		
		// Cutoff criterion based on Wald CI
		CI_lb = pvals[i_set] - Qz * sqrt(pvals[i_set] * (1 - pvals[i_set]) / n0_perm);
		if ( CI_lb <= (cutoff + tol)){
			// For p-value small enough, increase precision by adding more Y reps
			// For the remaining replicates (excluding obs Y)
			for( i_y = (n0_perm + 2); i_y <= (n_perm + 1); i_y++){
				
				// Step 1: For each variant, compute weight w_i
				for (i = 1; i <= n_marker; i++){
					for (m_u = 0, n_u = 0, fam = 1; fam <= n_fam; fam++){
						genos = family[fam].genos;

						for (j = 1; j <= n_person; j++){
							if (family[fam].trait[i_y][j] < tol){
								founder_allele = genos[j][0] + (n_f * 2)*(fam - 1);
								m_u += (haplo[line_in_haplo[i_set][founder_allele]][marknum[i_set][i]] == 1) ? 1 : 0;
								founder_allele = genos[j][1] + (n_f * 2)*(fam - 1);
								m_u += (haplo[line_in_haplo[i_set][founder_allele]][marknum[i_set][i]] == 1) ? 1 : 0;
								n_u++;
							}
						}
					}
					w[i] = (1.0 + m_u) / (2.0 * (1 + n_u));
					w[i] = sqrt(n_total * w[i] * (1.0 - w[i]));
				}
				
				// Step 2: For each individual, compute score
				for (k = 1, fam = 1; fam <= n_fam; fam++){
					genos = family[fam].genos;

					for (j = 1; j <= n_person; j++, k++){
						for (x = 0, i = 1; i <= n_marker; i++){
							founder_allele = genos[j][0] + (n_f * 2)*(fam - 1);
							Gi = (haplo[line_in_haplo[i_set][founder_allele]][marknum[i_set][i]] == 1) ? 1 : 0;
							founder_allele = genos[j][1] + (n_f * 2)*(fam - 1);
							Gi += (haplo[line_in_haplo[i_set][founder_allele]][marknum[i_set][i]] == 1) ? 1 : 0;
							Gi /= w[i];
							x += Gi;
						}
						score[k] = x, score_ind[k] = k;
					}
				}

				// Step 3: Get the rank for all individuals and compute the sum of the ranks for cases
				getRank(score, score_ind, ranks, n_total);
				for (x = 0, k = 1, fam = 1; fam <= n_fam; fam++){
					for (j = 1; j <= n_person; j++, k++){
						x += ranks[k] * family[fam].trait[i_y][j];
					}
				}
				// Compute the p-value from the new rep by comparing zVT between Y and Y rep
				if ((zWS + tol) < x) p_vals++;
				else if (fabs(zWS - x) <= tol ) p_vals += (rand1(seed) < .5) ? 1 : 0;
			}
			pvals[i_set] = (p_vals + 1.0) / (n_perm + 1.0);
		}
	}
	free_dvector(w, 1, n_marker);
	free_dvector(score, 1, n_total);
	free_ivector(score_ind, 1, n_total);
	free_dvector(ranks, 1, n_total);
}


int cmp(const void *a, const void *b){
	double aa = *(double*)a, bb = *(double*)b;
	return (aa > bb) - (aa < bb);
}

void getRank(double *A, int *B, double *rank, int n_a){
	int i = 1, j, k;
	double sum_rank;
	
	array = A;
	qsort(&B[1], n_a, sizeof(B[1]), cmp_rank);
	
	while(i < n_a){
		rank[ B[i] ] = i;  // Index i
		sum_rank = i, k = 1; 
		while( (i < n_a) && (A[ B[i] ] == A[ B[i + 1] ]) ) sum_rank += ++i, k++; 
		i++;
		if (k > 1){
			for(j = 1; j <= k; j++) rank[ B[i-j] ] = sum_rank / k ; // For ties, take average of ranks
		}
		if ((i == n_a) && (A[ B[i] ] != A[ B[i - 1] ])) rank[ B[i] ] = i; // For the last entry if it's unique
	}
}

int cmp_rank(const void *a, const void *b){
	int aa = *(int*)a, bb = *(int*)b;
	return (array[aa]> array[bb]) - (array[aa] < array[bb]);
}
