/*

    _dtram.h - dTRAM implementation in C (header file)

    author: Christoph Wehmeyer <christoph.wehmeyer@fu-berlin.de>

*/

#ifndef PYTRAM_DTRAM
#define PYTRAM_DTRAM

#include <math.h>

void _nu_K_ij_equation(
    double *nu_K_i,
    double *gamma_K_i,
    double *pi_i,
    int *C_K_ij,
    int n_therm_states,
    int n_markov_states,
    double *new_nu_K_i
);

void _pi_i_equation(
    double *nu_K_i,
    double *gamma_K_i,
    double *pi_i,
    int *C_K_ij,
    int n_therm_states,
    int n_markov_states,
    double *new_pi_i
);

void _p_K_ij_equation(
    double *nu_K_i,
    double *gamma_K_i,
    double *pi_i,
    int *C_K_ij,
    int n_therm_states,
    int n_markov_states,
    double *p_K_ij
);

#endif
