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

def dtram( C_K_ij, gamma_K_ij, maxiter=100, ftol=1.0E-5, verbose=False ):
    r"""

    """
    # try to create the DTRAM object
    try:
        dtram_obj = DTRAM( C_K_ij, gamma_K_ij )
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


