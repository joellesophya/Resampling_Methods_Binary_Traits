#ifndef PERMUTE_H
#define PERMUTE_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "simulation.h"
#include "cholesky.h"
#include "svdcomp.h"
#include "nrutil.h"

void permute_pre(struct FAMILY *family, struct DATA_STRUCT data_struct, struct PHI_SVD phi_svd, struct TMP_PERM tmp_perm;);
void permute(struct FAMILY *family, struct DATA_STRUCT data_struct, struct TMP_PERM tmp_perm);
int cmpfunc (const void * a, const void * b);

#endif
