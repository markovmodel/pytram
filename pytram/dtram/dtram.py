r"""

======================
dTRAM estimator module
======================

.. moduleauthor:: Christoph Wehmeyer <christoph.wehmeyer@fu-berlin.de>

"""

import numpy as np
from ..estimator import Estimator, NotConvergedWarning, ExpressionError
from .ext import nu_K_i_equation_tagged, pi_i_equation_tagged, p_K_ij_equation_tagged
from .ext import log_nu_K_i_setter_lse, log_nu_K_i_equation_lse, f_i_equation_lse, p_K_ij_equation_lse



####################################################################################################
#
#   DTRAM ESTIMATOR CLASS
#
####################################################################################################

class DTRAM( Estimator ):
    r"""
    The dTRAM estimator class
    """
    def __init__( self, C_K_ij, b_K_i, tagged=True ):
        r"""
        Initialize the DTRAM object
        
        Parameters
        ----------
        C_K_ij : numpy.ndarray( shape=(T,M,M), dtype=numpy.intc )
            transition counts between the M discrete Markov states for each of the T thermodynamic ensembles
        b_K_i : numpy.ndarray( shape=(T,M), dtype=numpy.float64 )
            bias energies in the T thermodynamic and M discrete Markov states
        """
        super( DTRAM, self ).__init__( C_K_ij )
        self.tagged = tagged
        if tagged:
            # this check raises an exception if b_K_i is not usable
            if self._check_b_K_i( b_K_i ):
                self._gamma_K_i = np.exp( b_K_i.min() - b_K_i )
            # hard-coded initial guess for pi_i and nu_K_i
            self._pi_i = np.ones( shape=(self.n_markov_states,), dtype=np.float64 ) / float( self.n_markov_states )
            self._nu_K_i = C_K_ij.sum( axis=2 ).astype( np.float64 )
            # 'private' storage variable initialization
            self._f_K = None
            self._pi_K_i = None
        else:
            pass
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
    def pi_i( self ):
        if self.tagged:
            return self._pi_i
        else:
            return None

    @property
    def pi_K_i( self ):
        if self.tagged:
            if self._pi_K_i is None:
                self._pi_K_i = self.f_K[:,np.newaxis] * self.pi_i[np.newaxis,:] * self.gamma_K_i
            return self._pi_K_i
        else:
            return None

    @property
    def f_K( self ):
        if self.tagged:
            if self._f_K is None:
                self._f_K = 1.0 / np.dot( self.gamma_K_i, self.pi_i )
            return self._f_K
        else:
            return None

    @property
    def nu_K_i( self ):
        if self.tagged:
            return self._nu_K_i
        else:
            return None

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
        maxiter : int
            maximum number of self-consistent-iteration steps
        ftol : float (> 0.0)
            convergence criterion based on the max relative change in an self-consistent-iteration step
        verbose : boolean
            writes convergence information to stdout during the self-consistent-iteration cycle
        """
        if self.tagged:
            # reset internal storage variables
            self._f_K = None
            self._pi_K_i = None
            finc = None
            if verbose:
                print "# %25s %25s" % ( "[iteration step]", "[relative increment]" )
            # start the iteration loop
            for i in xrange( maxiter ):
                # iterate nu_K_i
                tmp_nu_K_i = np.copy( self.nu_K_i )
                nu_K_i_equation_tagged( tmp_nu_K_i, self.gamma_K_i, self.pi_i, self.C_K_ij, self._nu_K_i )
                # iterate pi_i
                tmp_pi_i = np.copy( self.pi_i )
                pi_i_equation_tagged( self.nu_K_i, self.gamma_K_i, tmp_pi_i, self.C_K_ij, self._pi_i )
                # normalize pi_i
                self._pi_i /= self.pi_i.sum()
                # compute the relative change of pi_i
                nonzero = tmp_pi_i.nonzero()
                div = np.ones( shape=(self.n_markov_states,), dtype=np.float64 )
                div[nonzero] = tmp_pi_i[nonzero]
                finc = np.max( np.abs( tmp_pi_i - self.pi_i ) / div )
                # write out progress if requested
                if verbose:
                    print " %25d %25.12e" % ( i+1, finc )
                # break loop if we're converged
                if finc < ftol:
                    break
            # complain if we're not yet converged
            if finc > ftol:
                raise NotConvergedWarning( "DTRAM", finc )
        else:
            pass

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
        if self.tagged:
            p_K_ij_equation_tagged( self.nu_K_i, self.gamma_K_i, self.pi_i, self.C_K_ij, p_K_ij )
        else:
            pass
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
    #   gamma_K_i sanity checks and getter
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
            raise ExpressionError( "b_K_i", "unmatching number of thermodynamic states (%d,%d)" % (b_K_i.shape[0], self.n_therm_states) )
        if b_K_i.shape[1] != self.n_markov_states:
            raise ExpressionError( "b_K_i", "unmatching number of markov states (%d,%d)" % (b_K_i.shape[1], self.n_markov_states) )
        if np.float64 != b_K_i.dtype:
            raise ExpressionError( "b_K_i", "invalid dtype (%s)" % str( b_K_i.dtype ) )
        return True

    @property
    def gamma_K_i( self ):
        if self.tagged:
            return self._gamma_K_i
        else:
            return None

