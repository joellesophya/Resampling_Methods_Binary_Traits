#ifndef READ_H
#define READ_H

#include "simulation.h"

void readsimparam(char *filename, struct SIM_PARAM *param);
void pednum(char* filename, struct FAMILY *family, int n_cov, int n_marker);
void readped(char* filename, struct FAMILY *family);
void readkininb(char *filename, struct FAMILY *family);
void readhaplo(char *filename, int **haplotypes, int n_sites);


#endif
