#ifndef ASSOTEST_H
#define ASSOTEST_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nrutil.h"
#include "cholesky.h"
#include "simulation.h"

void wst(struct FAMILY *family, struct RV_INFO rv, int i_set, struct DATA_STRUCT data_struct, struct STATS stats, int perm);
void VT(struct FAMILY *family, struct RV_INFO rv, int i_set, struct DATA_STRUCT data_struct, struct PHI_SVD phi_svd, struct STATS stats, int perm);
int cmp(const void *a, const void *b);

#endif
