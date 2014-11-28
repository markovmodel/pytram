r"""

========================
API for the TRAM package
========================

.. moduleauthor:: Christoph Wehmeyer <christoph.wehmeyer@fu-berlin.de>
"""

from .dtram import DTRAM
from . import NotConvergedWarning, ExpressionError



####################################################################################################
#
#   dTRAM API function
#
####################################################################################################

def dtram( C_K_ij, gamma_K_i, maxiter=100, ftol=1.0E-5, verbose=False ):
    r"""
    The dTRAM API function
    
    Parameters
    ----------
    C_K_ij : numpy.ndarray( shape=(T,M,M), dtype=numpy.intc )
        transition counts between the M discrete Markov states for each of the T thermodynamic ensembles
    gamma_K_i : numpy.ndarray( shape=(T,M), dtype=numpy.float64 )
        conversion factors between the T thermodynamic and M discrete Markov states
    maxiter : int
        maximum number of SCF iteration steps during the optimisation of the stationary probabilities
    ftol : float (> 0.0)
        convergence criterion based on the max relative change in an SCF iteration step
    verbose : boolean
        writes convergence information to stdout during the SCF cycle
    
    Returns
    -------
    dtram_obj : object
        dTRAM estimator object with optimised unbiased stationary probabilities
    """
    # try to create the DTRAM object
    try:
        dtram_obj = DTRAM( C_K_ij, gamma_K_i )
    except ExpressionError, e:
        print "# ERROR ############################################################################"
        print "# Your input was faulty!"
        print "# The < %s > object is malformed: %s" % ( e.expression, e.msg )
        print "# ABORTING #########################################################################"
        return None
    # try to converge the stationary probabilities
    try:
        dtram_obj.scf_iteration( maxiter=maxiter, ftol=ftol, verbose=verbose )
    except NotConvergedWarning, e:
        print "# WARNING ##########################################################################"
        print "# dTRAM did not converge within %d steps!" % maxiter
        print "# The last relative increment was %.6e." % e.relative_increment
        print "# You should run the SCF iteration method again."
        print "# USE RESULTS WITH CARE ############################################################"
    finally:
        return dtram_obj


