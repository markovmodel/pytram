/*

    _dtram.c - dTRAM implementation in C

    author: Christoph Wehmeyer <christoph.wehmeyer@fu-berlin.de>

*/

#include "_dtram.h"



void _log_nu_K_i_setter_lse(
    double *log_nu_K_i,
    int *C_K_ij,
    int n_therm_states,
    int n_markov_states
)
{
    int i, j, K;
    int MM = n_markov_states * n_markov_states, KMM;
    int sum;
    for( K=0; K<n_therm_states; ++K )
    {
        KMM = K * MM;
        for( i=0; i<n_markov_states; ++i )
        {
            sum = 0;
            for( j=0; j<n_markov_states; ++j )
                sum += C_K_ij[ KMM + i*n_markov_states + j ];
            log_nu_K_i[ K * n_markov_states + i ] = log( PYTRAM_DTRAM_PRIOR + (double) sum );
        }
    }
}

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
{
    int i, j, K;
    int MM = n_markov_states * n_markov_states, Ki, Kj;
    int CK, CKij, CKji;
    double divisor;
    for( i=0; i<n_markov_states; ++i )
    {
        for( K=0; K<n_therm_states; ++K )
        {
            Ki = K*n_markov_states+i;
            for( j=0; j<n_markov_states; ++j )
            {
                CKij = C_K_ij[K*MM+i*n_markov_states+j];
                CKji = C_K_ij[K*MM+j*n_markov_states+i];
                /* special case: most variables cancel out, here */
                if( i == j )
                {
                    scratch_j[j] = ( 0 == CKij ) ? PYTRAM_DTRAM_LOG_PRIOR : log( PYTRAM_DTRAM_PRIOR + (double) CKij );
                    continue;
                }
                CK = CKij + CKji;
                Kj = K*n_markov_states+j;
                /* special case */
                if( 0 == CK )
                {
                    if( ( -INFINITY == log_nu_K_i[Ki] ) && ( -INFINITY == log_nu_K_i[Kj] ) )
                    {
                        scratch_j[j] = -_logsumexp_pair( 0.0, f_i[j] - f_i[i] + b_K_i[Kj] - b_K_i[Ki] );
                        printf( "####### WARNING ####### scratch_j=%f\n", scratch_j[j] );
                    }
                    else
                        scratch_j[j] = -INFINITY;
                    continue;
                }
                /* regular case */
                divisor = _logsumexp_pair(
                        log_nu_K_i[Kj] - f_i[i] - b_K_i[Ki],
                        log_nu_K_i[Ki] - f_i[j] - b_K_i[Kj]
                    );
                scratch_j[j] = log( (double) CK ) - b_K_i[Kj] - f_i[j] + log_nu_K_i[Ki] - divisor;
            }
            new_log_nu_K_i[Ki] = _logsumexp( scratch_j, n_markov_states );
        }
    }
}

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
{
    int i, j, K;
    int MM = n_markov_states * n_markov_states, KM = n_therm_states * n_markov_states, Ki, Kj;
    int CK, CKij, CKji, Ci;
    double divisor, norm;
    for( i=0; i<n_markov_states; ++i )
    {
        Ci = 0;
        for( K=0; K<n_therm_states; ++K )
        {
            Ki = K * n_markov_states + i;
            for( j=0; j<n_markov_states; ++j )
            {
                Kj = K * n_markov_states + j;
                CKij = C_K_ij[K*MM+i*n_markov_states+j];
                CKji = C_K_ij[K*MM+j*n_markov_states+i];
                /* add counts to Ci */
                Ci += CKji;
                /* special case: most variables cancel out, here */
                if( i == j )
                {
                    scratch_K_j[Kj] = ( 0 == CKij ) ? PYTRAM_DTRAM_LOG_PRIOR : log( PYTRAM_DTRAM_PRIOR + (double) CKij );
                    scratch_K_j[Kj] += f_i[i];
                    continue;
                }
                CK = CKij + CKji;
                /* special case */
                if( 0 == CK )
                {
                    scratch_K_j[Kj] = -INFINITY;
                    continue;
                }
                /* special case */ /* NaNs possible! CHECK THIS */
                /*if( -INFINITY == log_nu_K_i[Ki] )
                    continue;*/
                /* regular case */
                divisor = _logsumexp_pair(
                        log_nu_K_i[Kj] - f_i[i] - b_K_i[Ki],
                        log_nu_K_i[Ki] - f_i[j] - b_K_i[Kj]
                    );
                scratch_K_j[Kj] = log( (double) CK ) - b_K_i[Ki] + log_nu_K_i[Kj] - divisor;
            }
        }
        /* patch Ci and the total divisor together */
        new_f_i[i] = _logsumexp( scratch_K_j, KM ) - log( n_therm_states * PYTRAM_DTRAM_PRIOR + (double) Ci );
        scratch_j[i] = -new_f_i[i];
    }
    norm = _logsumexp( scratch_j, n_markov_states );
    for( i=0; i<n_markov_states; ++i )
        new_f_i[i] += norm;
}

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
{
    int i, j, K;
    int MM = n_markov_states * n_markov_states, KMM, Ki, Kj, ij, ji;
    int CK;
    double divisor, pKi;
    for( K=0; K<n_therm_states; ++K )
    {
        KMM = K*MM;
        for( i=0; i<n_markov_states; ++i )
        {
            pKi = 0.0;
            Ki = K*n_markov_states+i;
            for( j=0; j<n_markov_states; ++j )
            {
                /* special case: we compute the diagonal elements later */
                if( i == j )
                {
                    scratch_j[j] = -INFINITY;
                    continue;
                }
                ij = i*n_markov_states+j;
                ji = j*n_markov_states+i;
                p_K_ij[KMM+ij] = 0.0;
                CK = C_K_ij[KMM+ij] + C_K_ij[KMM+ji];
                /* special case: this element is zero */
                if( 0 == CK )
                {
                    scratch_j[j] = -INFINITY;
                    continue;
                }
                /* regular case */
                Kj = K*n_markov_states+j;
                divisor = _logsumexp_pair(
                        log_nu_K_i[Kj] - f_i[i] - b_K_i[Ki],
                        log_nu_K_i[Ki] - f_i[j] - b_K_i[Kj]
                    );
                scratch_j[j] =  log( (double) CK ) - f_i[j] - b_K_i[Kj] - divisor;
                p_K_ij[KMM+ij] = exp( scratch_j[j] );
            }
            /* compute the diagonal elements from the other elements in this line */
            p_K_ij[KMM+i*n_markov_states+i] = 1.0 - exp( _logsumexp( scratch_j, n_markov_states ) );
        }
    }
}








/* *************************************** TAGGED VERSION *************************************** */



void _nu_K_i_equation_tagged(
    double *nu_K_i,
    double *gamma_K_i,
    double *pi_i,
    int *C_K_ij,
    int n_therm_states,
    int n_markov_states,
    double *new_nu_K_i
)
{
    int i, j, K;
    int MM = n_markov_states * n_markov_states, Ki, Kj;
    int CK, CKij, CKji;
    double divisor;
    for( i=0; i<n_markov_states; ++i )
    {
        for( K=0; K<n_therm_states; ++K )
        {
            Ki = K*n_markov_states+i;
            /* initialize new_nu_K_i at zero */
            new_nu_K_i[Ki] = 0.0;
            for( j=0; j<n_markov_states; ++j )
            {
                CKij = C_K_ij[K*MM+i*n_markov_states+j];
                CKji = C_K_ij[K*MM+j*n_markov_states+i];
                /* special case: most variables cancel out, here */
                if( i == j )
                {
                    new_nu_K_i[Ki] += ( 0 == CKij ) ? 0.0 : (double) CKij;
                    /* new_nu_K_i[Ki] += 1.0E-10 + (double) CKij; */
                    continue;
                }
                CK = CKij + CKji;
                /* special case: we can skip this j */
                if( 0 == CK )
                    continue;
                /* special case: we can skip this j */
                if( 0.0 == nu_K_i[Ki] )
                    continue;
                /* regular case */
                Kj = K*n_markov_states+j;
                divisor = gamma_K_i[Kj] * pi_i[j] * nu_K_i[Ki];
                divisor += ( 0.0 == nu_K_i[Kj] ) ? 0.0 : gamma_K_i[Ki] * pi_i[i] * nu_K_i[Kj];
                new_nu_K_i[Ki] += nu_K_i[Ki] * (double) CK * gamma_K_i[Kj] * pi_i[j] / divisor;
            }
        }
    }
}

void _pi_i_equation_tagged(
    double *nu_K_i,
    double *gamma_K_i,
    double *pi_i,
    int *C_K_ij,
    int n_therm_states,
    int n_markov_states,
    double *new_pi_i
)
{
    int i, j, K;
    int MM = n_markov_states * n_markov_states, Ki, Kj;
    int CK, CKij, CKji, Ci;
    double divisor, div_i;
    for( i=0; i<n_markov_states; ++i )
    {
        Ci = 0;
        div_i = 0.0;
        for( K=0; K<n_therm_states; ++K )
        {
            Ki = K * n_markov_states + i;
            for( j=0; j<n_markov_states; ++j )
            {
                CKij = C_K_ij[K*MM+i*n_markov_states+j];
                CKji = C_K_ij[K*MM+j*n_markov_states+i];
                /* add counts to Ci */
                Ci += CKji;
                /* special case: most variables cancel out, here */
                if( i == j )
                {
                    div_i += ( 0 == CKij ) ? 0.0 : (double) CKij / pi_i[i];
                    /* new_pi_i[i] += ( 1.0E-10 + (double) CKij ) / pi_i[i]; */
                    continue;
                }
                CK = CKij + CKji;
                /* special case: we can skip this j */
                if( 0 == CK )
                    continue;
                Kj = K * n_markov_states + j;
                /* special case: we can skip this j */
                if( 0.0 == nu_K_i[Kj] )
                    continue;
                /* regular case */
                divisor = gamma_K_i[Ki] * pi_i[i] * nu_K_i[Kj];
                divisor += ( 0.0 == nu_K_i[Ki]) ? 0.0 : gamma_K_i[Kj] * pi_i[j] * nu_K_i[Ki];
                div_i += (double) CK * gamma_K_i[Ki] * nu_K_i[Kj] / divisor;
            }
        }
        /* patch Ci and the total divisor together */
        new_pi_i[i] = (double) Ci / div_i;
    }
}

void _p_K_ij_equation_tagged(
    double *nu_K_i,
    double *gamma_K_i,
    double *pi_i,
    int *C_K_ij,
    int n_therm_states,
    int n_markov_states,
    double *p_K_ij
)
{
    int i, j, K;
    int MM = n_markov_states * n_markov_states, KMM, Ki, Kj, ij, ji;
    int CK;
    double divisor, pKi;
    for( K=0; K<n_therm_states; ++K )
    {
        KMM = K*MM;
        for( i=0; i<n_markov_states; ++i )
        {
            pKi = 0.0;
            Ki = K*n_markov_states+i;
            for( j=0; j<n_markov_states; ++j )
            {
                /* special case: we compute the diagonal elements later */
                if( i == j )
                    continue;
                ij = i*n_markov_states+j;
                ji = j*n_markov_states+i;
                p_K_ij[KMM+ij] = 0.0;
                CK = C_K_ij[KMM+ij] + C_K_ij[KMM+ji];
                /* special case: this element is zero */
                if( 0 == CK )
                    continue;
                /* regular case */
                Kj = K*n_markov_states+j;
                divisor  = ( 0.0 == nu_K_i[Kj] ) ? 0.0 : gamma_K_i[Ki] * pi_i[i] * nu_K_i[Kj];
                divisor += ( 0.0 == nu_K_i[Ki] ) ? 0.0 : gamma_K_i[Kj] * pi_i[j] * nu_K_i[Ki];
                p_K_ij[KMM+ij] = CK * gamma_K_i[Kj] * pi_i[j] / divisor;
                pKi += p_K_ij[KMM+ij];
            }
            /* compute the diagonal elements from the other elements in this line */
            p_K_ij[KMM+i*n_markov_states+i] = 1.0 - pKi;
        }
    }
}




