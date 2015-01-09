r"""

========================
API for the TRAM package
========================

.. moduleauthor:: Christoph Wehmeyer <christoph.wehmeyer@fu-berlin.de>
.. moduleauthor:: Antonia Mey <antonia.mey@fu-berlin.de>
"""

from .dtram import DTRAM
from .xtram import XTRAM
from . import NotConvergedWarning, ExpressionError



####################################################################################################
#
#   dTRAM API function using the mathematical expressions at input
#
####################################################################################################

def dtram_me( C_K_ij, b_K_i, maxiter=100, ftol=1.0E-5, verbose=False ):
    r"""
    The dTRAM API function
    
    Parameters
    ----------
    C_K_ij : numpy.ndarray( shape=(T,M,M), dtype=numpy.intc )
        transition counts between the M discrete Markov states for each of the T thermodynamic ensembles
    b_K_i : numpy.ndarray( shape=(T,M), dtype=numpy.float64 )
        reduced bias energies at the T thermodynamic and M discrete Markov states
    maxiter : int
        maximum number of SCF iteration steps during the optimisation of the stationary probabilities
    ftol : float (> 0.0)
        convergence criterion based on the max relative change in an self-consistent-iteration step
    verbose : boolean
        writes convergence information to stdout during the self-consistent-iteration cycle
    
    Returns
    -------
    dtram_obj : object
        dTRAM estimator object with optimised unbiased stationary probabilities
    """
    # try to create the DTRAM object
    try:
        dtram_obj = DTRAM( C_K_ij, b_K_i )
    except ExpressionError, e:
        print "# ERROR ############################################################################"
        print "# Your input was faulty!"
        print "# The < %s > object is malformed: %s" % ( e.expression, e.msg )
        print "# ABORTING #########################################################################"
        raise
    # try to converge the stationary probabilities
    try:
        dtram_obj.sc_iteration( maxiter=maxiter, ftol=ftol, verbose=verbose )
    except NotConvergedWarning, e:
        print "# WARNING ##########################################################################"
        print "# dTRAM did not converge within %d steps!" % maxiter
        print "# The last relative increment was %.6e." % e.relative_increment
        print "# You should run the < sc_iteration > method again."
        print "# USE RESULTS WITH CARE ############################################################"
    finally:
        return dtram_obj



####################################################################################################
#
#   dTRAM API function from TRAMData
#
####################################################################################################

def dtram( tramdata, lag=1, sliding_window=True, maxiter=100, ftol=1.0E-5, verbose=False ):
    r"""
    The dTRAM API function
    
    Parameters
    ----------
    tramdata : object
        container/converter for TRAM input data
    lag : int
        specify the lag time for C_K_ij calculation
    sliding_window : boolean
        use sliding windows to calculate C_K_ij
    maxiter : int
        maximum number of SCF iteration steps during the optimisation of the stationary probabilities
    ftol : float (> 0.0)
        convergence criterion based on the max relative change in an self-consistent-iteration step
    verbose : boolean
        writes convergence information to stdout during the self-consistent-iteration cycle
    
    Returns
    -------
    dtram_obj : object
        dTRAM estimator object with optimised unbiased stationary probabilities
    """
    return dtram_me( tramdata.get_C_K_ij( lag ), tramdata.b_K_i, maxiter=maxiter, ftol=ftol, verbose=verbose )


####################################################################################################
#
#   xTRAM API function using expressions at input
#
####################################################################################################

def xtram_me( C_K_ij, u_I_x, T_x, M_x, N_K_i, N_K, maxiter=100, ftol=1.0E-5, verbose=False ):
    r"""
    The xTRAM API function
    
    Parameters
    ----------
    C_K_ij : numpy.ndarray( shape=(T,M,M), dtype=numpy.intc )
        transition counts between the M discrete Markov states for each of the T thermodynamic ensembles
    u_I_x : numpy.ndarray( shape=(T,L), dtype=numpy.float64 )
        biased energy evaluated at thermodynamic stated T over all sampled data of length L
    T_x : numpy.ndarray( shape=(L), dtype=numpy.intc )
        thermodynamic states over all sampled data of length L
    M_x : numpy.ndarray( shape=(L), dytpe=numpy.intc )
        Markov states over all sampled data of length L
    N_K_i : numpy.ndarray( shape=(T,M), dtype=numpy.intc )
        number of times a markov state (M) is seen in a thermodynamic state (T)
    N_K : numpy.ndarray( shape=(T), dtype=numpy.intc )
        number of observed samples for each thermodynamic state T
    maxiter : int
        maximum number of self consistent iteration steps during the optimisation of the stationary probabilities
    ftol : float (> 0.0)
        convergence criterion based on the max relative change in an self-consistent-iteration step
    verbose : boolean
        writes convergence information to stdout during the self-consistent-iteration cycle
    
    Returns
    -------
    xtram_obj : object
        xTRAM estimator object with optimised unbiased stationary probabilities
    """
    # try to create the DTRAM object
    try:
        xtram_obj = XTRAM( C_K_ij, u_I_x, T_x, M_x, N_K_i, N_K )
    except ExpressionError, e:
        print "# ERROR ############################################################################"
        print "# Your input was faulty!"
        print "# The < %s > object is malformed: %s" % ( e.expression, e.msg )
        print "# ABORTING #########################################################################"
        raise
    # try to converge the stationary probabilities
    try:
        xtram_obj.sc_iteration( maxiter=maxiter, ftol=ftol, verbose=verbose )
    except NotConvergedWarning, e:
        print "# WARNING ##########################################################################"
        print "# xTRAM did not converge within %d steps!" % maxiter
        print "# The last relative increment was %.6e." % e.relative_increment
        print "# You should run the < sc_iteration > method again."
        print "# USE RESULTS WITH CARE ############################################################"
    finally:
        return xtram_obj

  
####################################################################################################
#
#   xTRAM API function from TRAMData
#
####################################################################################################

def xtram( tramdata, lag=1, sliding_window=True, maxiter=100, ftol=1.0E-5, verbose=False ):
    r"""
    The xTRAM API function
    
    Parameters
    ----------
    tramdata : object
        container/converter for TRAM input data
    lag : int
        specify the lag time for C_K_ij calculation
    sliding_window : boolean
        use sliding windows to calculate C_K_ij
    maxiter : int
        maximum number of SCF iteration steps during the optimisation of the stationary probabilities
    ftol : float (> 0.0)
        convergence criterion based on the max relative change in an self-consistent-iteration step
    verbose : boolean
        writes convergence information to stdout during the self-consistent-iteration cycle
    
    Returns
    -------
    xtram_obj : object
        xTRAM estimator object with optimised unbiased stationary probabilities
    """
    return xtram_me( tramdata.get_C_K_ij( lag ), tramdata.u_I_x, tramdata.T_x, tramdata.M_x, tramdata.N_K_i, tramdata.N_K, maxiter=maxiter, ftol=ftol, verbose=verbose )
