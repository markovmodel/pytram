r"""

======================
dTRAM estimator module
======================

.. moduleauthor:: Christoph Wehmeyer <christoph.wehmeyer@fu-berlin.de>

"""

import numpy as np
from ..estimator import Estimator, NotConvergedWarning, ExpressionError
from .ext import log_nu_K_i_setter, log_nu_K_i_equation, f_i_equation, p_K_ij_equation, f_K_equation



####################################################################################################
#
#   DTRAM ESTIMATOR CLASS
#
####################################################################################################

class DTRAM( Estimator ):
    r"""
    This is the DTRAM class

    Parameters
    ----------
    C_K_ij : numpy.ndarray( shape=(T,M,M), dtype=numpy.intc )
        transition counts between the M discrete Markov states for each of the T thermodynamic ensembles
    b_K_i : numpy.ndarray( shape=(T,M), dtype=numpy.float64 )
        bias energies in the T thermodynamic and M discrete Markov states
    """

    def __init__( self, C_K_ij, b_K_i ):
        super( DTRAM, self ).__init__( C_K_ij )
        # this check raises an exception if b_K_i is not usable
        if self._check_b_K_i( b_K_i ):
            self._b_K_i = b_K_i
        # hard-coded initial guess for pi_i and nu_K_i
        self._f_i = np.zeros( shape=(self.n_markov_states,), dtype=np.float64 )
        self._log_nu_K_i = np.zeros( shape=(self.n_therm_states,self.n_markov_states), dtype=np.float64 )
        log_nu_K_i_setter( self._log_nu_K_i, self.C_K_ij )
        # citation information
        self.citation = [
                "Statistically optimal analysis of state-discretized trajectory data",
                "from multiple thermodynamic states;",
                "Hao Wu, Antonia S.J.S. Mey, Edina Rosta, and Frank Noe",
                "J. Chem. Phys. 141, 214106 (2014)"
            ]

    ############################################################################
    #
    #   override getters for stationary properties
    #
    ############################################################################

    @property
    def f_i( self ):
        return self._f_i

    @property
    def f_K_i( self ):
        return self.f_K[:,np.newaxis] + self._b_K_i + self._f_i[np.newaxis,:]

    @property
    def f_K( self ):
        _f_K = np.zeros( shape=(self.n_therm_states,), dtype=np.float64 )
        scratch_j = np.zeros( shape=(self.n_markov_states,), dtype=np.float64 )
        f_K_equation( self._b_K_i, self._f_i, scratch_j, _f_K )
        return _f_K

    ############################################################################
    #
    #   getters for dTRAM-specific properties
    #
    ############################################################################

    @property
    def b_K_i( self ):
        return self._b_K_i

    @property
    def gamma_K_i( self ):
        return np.exp( -self.b_K_i )

    @property
    def nu_K_i( self ):
        return np.exp( self._log_nu_K_i )

    ############################################################################
    #
    #   self-consistent-iteration to converge pi_i
    #
    ############################################################################

    def sc_iteration( self, maxiter=100, ftol=1.0E-5, verbose=False ):
        r"""
        Run the self-consistent-iteration cycle to optimise the unbiased stationary probabilities (and Langrange multipliers)
        
        Parameters
        ----------
        maxiter : int (default=100)
            maximum number of self-consistent-iteration steps
        ftol : float (default=1.0E-5)
            convergence criterion based on the max relative change in an self-consistent-iteration step
        verbose : boolean (default=False)
            writes convergence information to stdout during the self-consistent-iteration cycle
        """
        scratch_K_j = np.zeros( shape=(self.n_therm_states,self.n_markov_states), dtype=np.float64 )
        scratch_j = np.zeros( shape=(self.n_markov_states,), dtype=np.float64 )
        if verbose:
            print "# %25s %25s" % ( "[iteration step]", "[increment]" )
        # start the iteration loop
        for i in xrange( maxiter ):
            # iterate log_nu_K_i
            tmp_log_nu_K_i = np.copy( self._log_nu_K_i )
            log_nu_K_i_equation( tmp_log_nu_K_i, self._b_K_i, self._f_i, self.C_K_ij, scratch_j, self._log_nu_K_i )
            # iterate f_i
            tmp_f_i = np.copy( self._f_i )
            f_i_equation( self._log_nu_K_i, self._b_K_i, tmp_f_i, self.C_K_ij, scratch_K_j, scratch_j, self._f_i )
            # compute the absolute change of f_i
            finc = np.max( np.abs( tmp_f_i - self._f_i ) )
            # write out progress if requested
            if verbose:
                print " %25d %25.12e" % ( i+1, finc )
            # break loop if we're converged
            if finc < ftol:
                break
        # complain if we're not yet converged
        if finc > ftol:
            raise NotConvergedWarning( "DTRAM", finc )

    ############################################################################
    #
    #   transition matrix estimation
    #
    ############################################################################

    def estimate_transition_matrices( self ):
        r"""
        Estimate the transition matrices for all thermodynamic states
        
        Returns
        -------
        p_K_ij : numpy.ndarray( shape=(T,M,M), dtype=numpy.float64 )
            the transition matrices for all thermodynamic states
        """
        p_K_ij = np.zeros( shape=self.C_K_ij.shape, dtype=np.float64 )
        scratch_j = np.zeros( shape=(self.n_markov_states,), dtype=np.float64 )
        p_K_ij_equation( self._log_nu_K_i, self._b_K_i, self._f_i, self.C_K_ij, scratch_j, p_K_ij )
        return p_K_ij

    def estimate_transition_matrix( self, I ):
        r"""
        Estimate the transition matrices for one thermodynamic state
        
        Parameters
        ----------
        I : int
            target thermodynamic state
        
        Returns
        -------
        p_K_ij[I] : numpy.ndarray( shape=(M,M), dtype=numpy.float64 )
            the transition matrix for the Ith thermodynamic state
        """
        return self.estimate_transition_matrices()[I,:,:]

    ############################################################################
    #
    #   gamma_K_i sanity checks
    #
    ############################################################################

    def _check_b_K_i( self, b_K_i ):
        if b_K_i is None:
            raise ExpressionError( "b_K_i", "is None" )
        if not isinstance( b_K_i, (np.ndarray,) ):
            raise ExpressionError( "b_K_i", "invalid type (%s)" % str( type( b_K_i ) ) )
        if 2 != b_K_i.ndim:
            raise ExpressionError( "b_K_i", "invalid number of dimensions (%d)" % b_K_i.ndim )
        if b_K_i.shape[0] != self.n_therm_states:
            raise ExpressionError( "b_K_i", "not matching number of thermodynamic states (%d,%d)" % (b_K_i.shape[0], self.n_therm_states) )
        if b_K_i.shape[1] != self.n_markov_states:
            raise ExpressionError( "b_K_i", "not matching number of markov states (%d,%d)" % (b_K_i.shape[1], self.n_markov_states) )
        if np.float64 != b_K_i.dtype:
            raise ExpressionError( "b_K_i", "invalid dtype (%s)" % str( b_K_i.dtype ) )
        return True

