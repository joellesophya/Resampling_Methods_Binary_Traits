#ifndef GENOSIM_H
#define GENOSIM_H

void drop(int nfou, int nnon, int **parentmat, int **geno,  long *seed);
void codealleles(int n_founder_allele, int *state, long *seed);
void recode(int n_person, int *state, int **geno_in, int **geno_out);
void run_geno_sim(struct FAMILY family,int n_marker, long *seed);

#endif
