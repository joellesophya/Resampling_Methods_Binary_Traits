#ifndef TRAIT_SIM_H
#define TRAIT_SIM_H
#include "simulation.h"

void run_trait_sim(struct FAMILY family, double *mles, int n_cov, int i_perm, long *seed);

void agesexeffect(double *mu, double **cov, int n_person, double *mles, int n_cov);
void mvn_ran_obs(double *mu, int n_person, double **phi, double add_var, long *seed);
void omega_fun(double **omega, int n_person, double **kin, double add_var);
void compute_trait_logistic(double *trait, double *pi, int n_person, long *seed);
// void compute_trait_liability(int *trait, double *L, int n_person, struct SIM_PARAM param, long *seed);

//int countallele(int *genotype);

#endif
