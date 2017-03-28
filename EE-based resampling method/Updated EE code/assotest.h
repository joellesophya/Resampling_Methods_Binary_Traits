#ifndef ASSOTEST_H
#define ASSOTEST_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nrutil.h"
#include "cholesky.h"
#include "simulation.h"

void ws(struct FAMILY *family, struct RV_INFO rv, int n_set, struct DATA_STRUCT data_struct, double *pvals, int n_perm, long *seed);
void VT(struct FAMILY *family, struct RV_INFO rv, int n_set, struct DATA_STRUCT data_struct, double *pvals, int n_perm, long *seed);
void convert_trait(struct FAMILY *family, struct DATA_STRUCT data_struct);
int cmp(const void *a, const void *b);
void getRank(double *A, int *B, double *rank, int n_a);
int cmp_rank(const void *a, const void *b);


#endif
