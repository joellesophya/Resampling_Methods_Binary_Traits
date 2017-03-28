#include "permute_sim.h"

double *vec;

double* permute(struct FAMILY *family, struct DATA_STRUCT data_struct, struct STATS stats, int n_perm, int bin_Y_perm) {
	int i, j, fam, i_perm, reject, reject_fam;

	int n_person = family[1].n_person; //Assume same number of indiv for all families
	int n_total = data_struct.n_total;
	int n_fam = data_struct.n_fam;
	int n_cov = data_struct.n_cov;
	int nonzero = n_total - n_cov;

	double *y_perm, *mu, fam_dev, devs, tol = data_struct.tol;
	double *delta = stats.delta;
	double **GCV = stats.GCV_mat;
	double *sumsY = data_struct.sumsY;
	
	int *order_dec = ivector(1, n_total);
	double *Yperm = dvector(1, n_total);
	double *naff_ypi = dvector(1, n_fam); // Records number of cases in each fam in Y replicate
	
	for(i_perm = 1; i_perm <= n_perm; i_perm++){
		do {
			// Permute entries of delta
			shuffle(delta, nonzero);

			// Get Y replicate
			for (sumsY[i_perm + 1] = sumsY[i_perm + 2 + n_perm] = 0, fam = 1; fam <= n_fam; fam++){
				y_perm = family[fam].trait[i_perm + 1];
				mu = family[fam].mu_hat;
				
				for (i = 1; i <= n_person; i++){
					for (y_perm[i] = mu[i], j = 1; j <= nonzero; j++) y_perm[i] += GCV[(fam-1) * n_person + i][j] * delta[j];
					if(bin_Y_perm == 1){
						order_dec[(fam - 1) * n_person + i] = (fam - 1) * n_person + i;
						Yperm[(fam - 1) * n_person + i] = y_perm[i];
					}
					else{
						sumsY[i_perm + 1] += y_perm[i];
						sumsY[i_perm + 2 + n_perm] += y_perm[i] * y_perm[i];
					}
				}
			}
			
			if(bin_Y_perm == 1){ // Convert Y replicate to binary
				// Obtain the ranks of the Ys (in decreasing order)
				vec = Yperm;
				qsort(&order_dec[1], n_total, sizeof(order_dec[1]), cmp_dec_rank);
				// If order_dec > n-sumY then set to 1 o.w. 0 & record number of cases
				for (fam = 1; fam <= n_fam; fam++){
					y_perm = family[fam].trait[i_perm + 1];
					for (naff_ypi[fam] = 0, i = 1; i <= n_person; i++){
						y_perm[i] = (order_dec[(fam - 1) * n_person + i] > (n_total - sumsY[1] + tol)) ? 1 : 0;
						naff_ypi[fam] += y_perm[i];
					}
				}
				reject_fam = 0; // Set counter to 0
				
				for (devs = 0, fam = 1; fam <= n_fam; fam++){
					fam_dev = fabs(naff_ypi[fam] - family[fam].naff);
					if(fam_dev > (.3 * n_person)){ // 30% max deviations within families
						reject_fam++;
						break;
					}
					devs += fam_dev;
				}
				
				if((devs > (0.1 * n_total)) || reject_fam != 0) reject = 1; // 10% max deviations all across families
				else{
					reject = 0;
					sumsY[i_perm + 1] = sumsY[1]; // Same ncase with obs. Y
					sumsY[i_perm + 2 + n_perm] = sumsY[n_perm + 2];
				}
			} 
			else reject = 0;
		} while((reject == 1) && (bin_Y_perm == 1));
	}
	free_ivector(order_dec, 1, n_person);
	free_dvector(Yperm, 1, n_total);
	free_dvector(naff_ypi, 1, n_fam);
	return sumsY;
}

int cmp_dec_rank(const void *a, const void *b){  // For decreasing ranks
	int aa = *(int*)a, bb = *(int*)b;
	return -(vec[aa]> vec[bb]) + (vec[aa] < vec[bb]);
}
