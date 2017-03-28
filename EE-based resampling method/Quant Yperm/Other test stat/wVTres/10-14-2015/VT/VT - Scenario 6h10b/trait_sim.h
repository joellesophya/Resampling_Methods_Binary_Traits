#ifndef TRAIT_SIM_H
#define TRAIT_SIM_H
#include "simulation.h"

void run_trait_sim(struct FAMILY family, struct SIM_PARAM param, int model, long *seed);

void agesexeffect(double *mu, double **cov, int n_person, double age_effect, double sex_effect);
void mvn_ran_obs(double *mu, int n_person, double **phi, double **delta_7, struct SIM_PARAM param, int model, long *seed);
void omega_fun(double **omega, int n_person, double **kin, double **delta_7, double add_var, double err_var, double dom_var);
void compute_trait_logistic(double *trait, double *pi, int n_person, long *seed);
void compute_trait_liability(double *trait, double *L, int n_person, struct SIM_PARAM param, long *seed);


//int countallele(int *genotype);

#endif
