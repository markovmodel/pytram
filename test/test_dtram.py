################################################################################
#
#   test_dtram.py - testing the dTRAM estimator
#
#   author: Christoph Wehmeyer <christoph.wehmeyer@fu-berlin.de>
#
################################################################################

from nose.tools import assert_raises, assert_true
from pytram import DTRAM, ExpressionError, NotConvergedWarning
import numpy as np

def test_expression_error_None():
    assert_raises( ExpressionError, DTRAM, np.ones( shape=(2,3,3), dtype=np.intc ), None )

def test_expression_error_int():
    assert_raises( ExpressionError, DTRAM, np.ones( shape=(2,3,3), dtype=np.intc ), 5 )

def test_expression_error_list():
    assert_raises( ExpressionError, DTRAM, np.ones( shape=(2,3,3), dtype=np.intc ), [1,2] )

def test_expression_error_dim():
    assert_raises( ExpressionError, DTRAM, np.ones( shape=(2,3,3), dtype=np.intc ), np.ones( shape=(2,2,2), dtype=np.float64 ) )

def test_expression_error_markov():
    assert_raises( ExpressionError, DTRAM, np.ones( shape=(2,3,3), dtype=np.intc ), np.ones( shape=(2,2), dtype=np.float64 ) )

def test_expression_error_therm():
    assert_raises( ExpressionError, DTRAM, np.ones( shape=(2,3,3), dtype=np.intc ), np.ones( shape=(1,3), dtype=np.float64 ) )

def test_expression_error_int16():
    assert_raises( ExpressionError, DTRAM, np.ones( shape=(2,3,3), dtype=np.intc ), np.ones( shape=(2,3), dtype=np.int16 ) )

def test_expression_error_float32():
    assert_raises( ExpressionError, DTRAM, np.ones( shape=(2,3,3), dtype=np.intc ), np.ones( shape=(2,3), dtype=np.float32 ) )

def test_three_state_model():
    C_K_ij = np.array( [[[2358,29,0],[29,0,32],[0,32,197518]],[[16818,16763,0],[16763,0,16510],[0,16510,16635]]], dtype=np.intc )
    b_K_i = np.array( [ [ 0.0, 0.0, 0.0 ], [ 4.0, 0.0, 8.0 ] ], dtype=np.float64 )
    dtram = DTRAM( C_K_ij, b_K_i )
    assert_raises( NotConvergedWarning, dtram.sc_iteration, maxiter=1, ftol=1.0E-80, verbose=False )
    dtram.sc_iteration( maxiter=200000, ftol=1.0E-15, verbose=False )
    pi = np.array( [1.82026887e-02,3.30458960e-04,9.81466852e-01], dtype=np.float64 )
    T = np.array( [[9.90504397e-01,9.49560284e-03,0.0],[5.23046803e-01,0.0,4.76953197e-01],[0.0,1.60589690e-04,9.99839410e-01]], dtype=np.float64 )
    assert_true( np.max( np.abs( dtram.pi_i - pi ) ) < 1.0E-8 )
    assert_true( np.max( np.abs( dtram.estimate_transition_matrix( 0 ) - T ) ) < 1.0E-8 )




