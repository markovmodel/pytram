r"""

==============================
Basic estimator class for TRAM
==============================

.. moduleauthor:: Christoph Wehmeyer <christoph.wehmeyer@fu-berlin.de>
"""

import numpy as np



####################################################################################################
#
#   BASIC TRAM ESTIMATOR CLASS
#
####################################################################################################

class Estimator( object ):

    r"""
    Estimator is the parent class for all estimators

    Parameters
    ----------
    C_K_ij : numpy.ndarray( shape=(T,M,M), dtype=numpy.intc )
        transition counts between the M discrete Markov states for each of the T thermodynamic ensembles
    """
    def __init__( self, C_K_ij ):
        # this check raises an exception if C_K_ij is not usable
        if self._check_C_K_ij( C_K_ij ):
            self._C_K_ij = C_K_ij
        # if we reach this point, C_K_ij is save
        self._n_therm_states = C_K_ij.shape[0]
        self._n_markov_states = C_K_ij.shape[1]
        self.citation = []

    def cite( self, pre="" ):
        r"""
        Parameters
        ----------
        pre : string (default="")
            prepend string for printing citation string
        """
        for line in self.citation:
            print "%s%s" % ( pre, line )

    ############################################################################
    #
    #   C_K_ij sanity checks and getter
    #
    ############################################################################

    def _check_C_K_ij( self, C_K_ij ):
        if C_K_ij is None:
            raise ExpressionError( "C_K_ij", "is None" )
        if not isinstance( C_K_ij, (np.ndarray,) ):
            raise ExpressionError( "C_K_ij", "invalid type (%s)" % str( type( C_K_ij ) ) )
        if 3 != C_K_ij.ndim:
            raise ExpressionError( "C_K_ij", "invalid number of dimensions (%d)" % C_K_ij.ndim )
        if C_K_ij.shape[1] != C_K_ij.shape[2]:
            raise ExpressionError( "C_K_ij", "unmatching number of markov states (%d,%d)" % (C_K_ij.shape[1], C_K_ij.shape[2]) )
        if np.intc != C_K_ij.dtype:
            raise ExpressionError( "C_K_ij", "invalid dtype (%s)" % str( C_K_ij.dtype ) )
        if not np.all( C_K_ij.sum( axis=(0,2) ) > 0 ):
            raise ExpressionError( "C_K_ij", "contains unvisited states" )
        # TODO: strong connectivity check?
        return True

    @property
    def C_K_ij( self ):
        return self._C_K_ij

    ############################################################################
    #
    #   compute the TRAM log likelihood
    #
    ############################################################################

    def estimate_log_L_TRAM( self, C_K_ij, p_K_ij ):
        nonzero = C_K_ij.nonzero()
        return np.sum( C_K_ij[nonzero] * np.log( p_K_ij[nonzero] ) )

    ############################################################################
    #
    #   getters for stationary properties
    #
    ############################################################################

    @property
    def n_therm_states( self ):
        return self._n_therm_states

    @property
    def n_markov_states( self ):
        return self._n_markov_states

    @property
    def pi_i( self ):
        return np.exp( -self.f_i )

    @property
    def pi_K_i( self ):
        return np.exp( -self.f_K_i )

    @property
    def f_K( self ):
        raise NotImplementedError( "Override in derived class!" )

    @property
    def f_K_i( self ):
        return -np.log( self.pi_K_i )

    @property
    def f_i( self ):
        return -np.log( self.pi_i )



####################################################################################################
#
#   ERROR CLASS FOR MALFORMED ESTIMATOR ARGUMENTS
#
####################################################################################################

class ExpressionError( Exception ):
    r"""
    Exception class for malformed expressions in the input
    """
    def __init__( self, expression, msg ):
        self.expression = expression
        self.msg = msg
    def __str__( self ):
        return "[%s] %s" % ( self.expression, self.msg )



####################################################################################################
#
#   WARNING CLASS FOR PREMATURELY TERMINATED SCF ITERATIONS
#
####################################################################################################

class NotConvergedWarning( Exception ):
    r"""
    Exception class for non-convergence of estimators
    """
    def __init__( self, estimator, increment ):
        self.estimator = estimator
        self.increment = increment
    def __str__( self ):
        return "[%s] only reached increment %.3e" % ( self.estimator, self.increment )
