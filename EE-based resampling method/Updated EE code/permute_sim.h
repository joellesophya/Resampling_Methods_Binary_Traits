#ifndef PERMUTE_H
#define PERMUTE_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "simulation.h"
#include "sim_pheno_geno.h"
#include "nrutil.h"

double* permute(struct FAMILY *family, struct DATA_STRUCT data_struct, struct STATS stats, int n_perm, int bin_Y_perm);
int cmp_dec_rank(const void *a, const void *b);

#endif
