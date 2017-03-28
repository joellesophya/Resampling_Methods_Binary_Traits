#ifndef ASSOTEST_H
#define ASSOTEST_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nrutil.h"
#include "simulation.h"

void wst(struct FAMILY *family, int **haplotypes, struct DATA_STRUCT data_struct, struct STATS *stats);
void VT(struct FAMILY *family, struct RV_INFO rv, int i_set, struct DATA_STRUCT data_struct, struct STATS *stats, double *MAF, double *sorted_MAF, int *unique_MAF, double *threshold, double *z, int perm);
int cmp(const void *a, const void *b);

#endif
