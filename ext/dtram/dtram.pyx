################################################################################
#
#   dtram.pyx - dTRAM implementation in C (cython wrapper)
#
#   author: Christoph Wehmeyer <christoph.wehmeyer@fu-berlin.de>
#
################################################################################

import numpy as np
cimport numpy as np

cdef extern from "_dtram.h":
    void _log_nu_K_i_setter_lse(
            double *log_nu_K_i,
            int *C_K_ij,
            int n_therm_states,
            int n_markov_states
        )
    void _log_nu_K_i_equation_lse(
            double *log_nu_K_i,
            double *b_K_i,
            double *f_i,
            int *C_K_ij,
            int n_therm_states,
            int n_markov_states,
            double *scratch_j,
            double *new_log_nu_K_i
        )
    void _f_i_equation_lse(
            double *log_nu_K_i,
            double *b_K_i,
            double *f_i,
            int *C_K_ij,
            int n_therm_states,
            int n_markov_states,
            double *scratch_K_j,
            double *scratch_j,
            double *new_f_i
        )
    void _p_K_ij_equation_lse(
            double *log_nu_K_i,
            double *b_K_i,
            double *f_i,
            int *C_K_ij,
            int n_therm_states,
            int n_markov_states,
            double *scratch_j,
            double *p_K_ij
        )
    void _nu_K_i_equation_tagged(
            double *nu_K_i,
            double *gamma_K_i,
            double* pi_i,
            int *C_K_ij,
            int n_therm_states,
            int n_markov_states,
            double *new_nu_K_i
        )
    void _pi_i_equation_tagged(
            double *nu_K_i,
            double *gamma_K_i,
            double* pi_i,
            int *C_K_ij,
            int n_therm_states,
            int n_markov_states,
            double *new_pi_i
        )
    void _p_K_ij_equation_tagged(
            double *nu_K_i,
            double *gamma_K_i,
            double *pi_i,
            int *C_K_ij,
            int n_therm_states,
            int n_markov_states,
            double *p_K_ij
        )

def log_nu_K_i_setter_lse(
        np.ndarray[double, ndim=2, mode="c"] log_nu_K_i not None,
        np.ndarray[int, ndim=3, mode="c"] C_K_ij not None
    ):
    _log_nu_K_i_setter_lse(
            <double*> np.PyArray_DATA( log_nu_K_i ),
            <int*> np.PyArray_DATA( C_K_ij ),
            log_nu_K_i.shape[0],
            log_nu_K_i.shape[1]
        )

def log_nu_K_i_equation_lse(
        np.ndarray[double, ndim=2, mode="c"] log_nu_K_i not None,
        np.ndarray[double, ndim=2, mode="c"] b_K_i not None,
        np.ndarray[double, ndim=1, mode="c"] f_i not None,
        np.ndarray[int, ndim=3, mode="c"] C_K_ij not None,
        np.ndarray[double, ndim=1, mode="c"] scratch_j not None,
        np.ndarray[double, ndim=2, mode="c"] new_log_nu_K_i not None
    ):
    _log_nu_K_i_equation_lse(
            <double*> np.PyArray_DATA( log_nu_K_i ),
            <double*> np.PyArray_DATA( b_K_i ),
            <double*> np.PyArray_DATA( f_i ),
            <int*> np.PyArray_DATA( C_K_ij ),
            log_nu_K_i.shape[0],
            log_nu_K_i.shape[1],
            <double*> np.PyArray_DATA( scratch_j ),
            <double*> np.PyArray_DATA( new_log_nu_K_i )
        )

def f_i_equation_lse(
        np.ndarray[double, ndim=2, mode="c"] log_nu_K_i not None,
        np.ndarray[double, ndim=2, mode="c"] b_K_i not None,
        np.ndarray[double, ndim=1, mode="c"] f_i not None,
        np.ndarray[int, ndim=3, mode="c"] C_K_ij not None,
        np.ndarray[double, ndim=2, mode="c"] scratch_K_j not None,
        np.ndarray[double, ndim=1, mode="c"] scratch_j not None,
        np.ndarray[double, ndim=1, mode="c"] new_f_i not None
    ):
    _f_i_equation_lse(
            <double*> np.PyArray_DATA( log_nu_K_i ),
            <double*> np.PyArray_DATA( b_K_i ),
            <double*> np.PyArray_DATA( f_i ),
            <int*> np.PyArray_DATA( C_K_ij ),
            log_nu_K_i.shape[0],
            log_nu_K_i.shape[1],
            <double*> np.PyArray_DATA( scratch_K_j ),
            <double*> np.PyArray_DATA( scratch_j ),
            <double*> np.PyArray_DATA( new_f_i )
        )

def p_K_ij_equation_lse(
        np.ndarray[double, ndim=2, mode="c"] log_nu_K_i not None,
        np.ndarray[double, ndim=2, mode="c"] b_K_i not None,
        np.ndarray[double, ndim=1, mode="c"] f_i not None,
        np.ndarray[int, ndim=3, mode="c"] C_K_ij not None,
        np.ndarray[double, ndim=1, mode="c"] scratch_j not None,
        np.ndarray[double, ndim=3, mode="c"] p_K_ij not None
    ):
    _p_K_ij_equation_lse(
            <double*> np.PyArray_DATA( log_nu_K_i ),
            <double*> np.PyArray_DATA( b_K_i ),
            <double*> np.PyArray_DATA( f_i ),
            <int*> np.PyArray_DATA( C_K_ij ),
            log_nu_K_i.shape[0],
            log_nu_K_i.shape[1],
            <double*> np.PyArray_DATA( scratch_j ),
            <double*> np.PyArray_DATA( p_K_ij )
        )



########################################## TAGGED VERSION ##########################################



def nu_K_i_equation_tagged(
        np.ndarray[double, ndim=2, mode="c"] nu_K_i not None,
        np.ndarray[double, ndim=2, mode="c"] gamma_K_i not None,
        np.ndarray[double, ndim=1, mode="c"] pi_i not None,
        np.ndarray[int, ndim=3, mode="c"] C_K_ij not None,
        np.ndarray[double, ndim=2, mode="c"] new_nu_K_i not None
    ):
    _nu_K_i_equation_tagged(
            <double*> np.PyArray_DATA( nu_K_i ),
            <double*> np.PyArray_DATA( gamma_K_i ),
            <double*> np.PyArray_DATA( pi_i ),
            <int*> np.PyArray_DATA( C_K_ij ),
            nu_K_i.shape[0],
            nu_K_i.shape[1],
            <double*> np.PyArray_DATA( new_nu_K_i )
        )

def pi_i_equation_tagged(
        np.ndarray[double, ndim=2, mode="c"] nu_K_i not None,
        np.ndarray[double, ndim=2, mode="c"] gamma_K_i not None,
        np.ndarray[double, ndim=1, mode="c"] pi_i not None,
        np.ndarray[int, ndim=3, mode="c"] C_K_ij not None,
        np.ndarray[double, ndim=1, mode="c"] new_pi_i not None
    ):
    _pi_i_equation_tagged(
            <double*> np.PyArray_DATA( nu_K_i ),
            <double*> np.PyArray_DATA( gamma_K_i ),
            <double*> np.PyArray_DATA( pi_i ),
            <int*> np.PyArray_DATA( C_K_ij ),
            nu_K_i.shape[0],
            nu_K_i.shape[1],
            <double*> np.PyArray_DATA( new_pi_i )
        )

def p_K_ij_equation_tagged(
        np.ndarray[double, ndim=2, mode="c"] nu_K_i not None,
        np.ndarray[double, ndim=2, mode="c"] gamma_K_i not None,
        np.ndarray[double, ndim=1, mode="c"] pi_i not None,
        np.ndarray[int, ndim=3, mode="c"] C_K_ij not None,
        np.ndarray[double, ndim=3, mode="c"] p_K_ij not None
    ):
    _p_K_ij_equation_tagged(
            <double*> np.PyArray_DATA( nu_K_i ),
            <double*> np.PyArray_DATA( gamma_K_i ),
            <double*> np.PyArray_DATA( pi_i ),
            <int*> np.PyArray_DATA( C_K_ij ),
            nu_K_i.shape[0],
            nu_K_i.shape[1],
            <double*> np.PyArray_DATA( p_K_ij )
        )
