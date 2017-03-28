#ifndef SIMPHENO_H
#define SIMPHENO_H

#define MAXLEN 1024
#define MISSVAL -9.0
#define NUMTOL 1e-6
#define ANALYSIS_STATUS 1

#include "simulation.h"

void sim_pheno_geno(struct FAMILY *family, struct RV_INFO rv, struct DATA_STRUCT data_struct,int n_sites, long *seed);
void store_haploline(struct RV_INFO rv, int n_founders, int n_sites, int n_marker);
void shuffle(double *array, int n);
static int rand_int(int n);

#endif
