/*

    _xtram.h - xTRAM implementation in C (header file)

    author: Antonia Mey <antonia.mey@fu-berlin.de>

*/

#ifndef PYTRAM_XTRAM
#define PYTRAM_XTRAM

#include <math.h>
#include <stdlib.h>
#include <stdio.h>


typedef struct
{
	int i;
	int j;
	double value;
}sparse_x;



void _b_i_IJ_equation(
	int T_length, 
	int n_therm_states, 
	int n_markov_states,
	int *T_x, 
	int *M_x,
	int *N,  
	double *f,  
	double *w,
	double *u,
	double *b_i_IJ);
	
double _iterate_x(
	long n_entries,
	long pi_length,
	long maxiter,
	double ftol,
	int *C_i,
	int *C_j,
	double *C_ij,
	double *C_ji,
	double *x_row,
	double *c_column,
	double *pi);

void update_x( 
	double *x_row, 
	sparse_x *x, 
	int *C_i, 
	int *C_j, 
	double *C_ij, 
	double *C_ji, 
	double *c_column, 
	int L);
void update_x_row(int L, sparse_x *x, double *x_row, int x_row_l);
void compute_pi(double *pi, double *x_row, int l_pi);
double converged(double *pi_old, double *pi_new, int l_pi);
void printpi(double *pi, int l);


#endif
