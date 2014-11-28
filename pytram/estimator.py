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
    I am the parent class for all estimators!
    """
    def __init__( self, C_K_ij ):
        self.C_K_ij = C_K_ij
        self.n_therm_states = C_K_ij.shape[0]
        self.n_markov_states = C_K_ij.shape[1]

    ############################################################################
    #
    #   C_K_ij sanity checks and getter/setter
    #
    ############################################################################

    def _check_C_K_ij( self, C_K_ij ):
        if None == C_K_ij:
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

    @C_K_ij.setter
    def C_K_ij( self, C_K_ij ):
        self._C_K_ij = None
        if self._check_C_K_ij( C_K_ij ):
            self._C_K_ij = C_K_ij

    ############################################################################
    #
    #   compute the TRAM log likelihood
    #
    ############################################################################

    def estimate_log_L_TRAM( self, C_K_ij, p_K_ij ):
        nonzero = C_K_ij.nonzero()
        return np.sum( C_K_ij[nonzero] * np.log( p_K_ij[nonzero] ) )



####################################################################################################
#
#   ERROR CLASS FOR MALFORMED ESTIMATOR ARGUMENTS
#
####################################################################################################

class ExpressionError( Exception ):
    r"""
    Exception class for malformed exressions in the input
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
    def __init__( self, estimator, relative_increment ):
        self.estimator = estimator
        self.relative_increment = relative_increment
    def __str__( self ):
        return "[%s] only reached relative increment %.3e" % ( self.estimator, self.relative_increment )
