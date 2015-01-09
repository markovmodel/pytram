################################################################################
#
#   test_xtram.py - testing the xTRAM estimator
#
#   author: Antonia Mey <antonia.mey@fu-berlin.de>
#
################################################################################

from nose.tools import assert_raises, assert_true
from pytram import XTRAM, ExpressionError, NotConvergedWarning
import numpy as np


def test_expression_error_None():
    assert_raises( ExpressionError, XTRAM, np.ones( shape =(2,3,3), dtype=np.intc), None, np.ones( shape =(10), dtype=np.intc), np.ones( shape =(10), dtype=np.intc),np.ones( shape =(2,3), dtype=np.intc), np.ones( shape =(2), dtype=np.intc) )
    assert_raises( ExpressionError, XTRAM, np.ones( shape =(2,3,3), dtype=np.intc), np.ones( shape =(2,10), dtype=np.float64 ),None, np.ones( shape =(10), dtype=np.intc), np.ones( shape =(2,3), dtype=np.intc), np.ones( shape =(2), dtype=np.intc) )
    assert_raises( ExpressionError, XTRAM, np.ones( shape =(2,3,3), dtype=np.intc), np.ones( shape =(2,10), dtype=np.float64), np.ones( shape =(10), dtype=np.intc),None, np.ones( shape =(2,3), dtype=np.intc), np.ones( shape =(2), dtype=np.intc) )
    assert_raises( ExpressionError, XTRAM, np.ones( shape =(2,3,3), dtype=np.intc), np.ones( shape =(2,10), dtype=np.float64), np.ones( shape =(10), dtype=np.intc),np.ones( shape =(10), dtype=np.intc), None, np.ones( shape =(), dtype=np.intc) )
    assert_raises( ExpressionError, XTRAM, np.ones( shape =(2,3,3), dtype=np.intc), np.ones( shape =(2,10), dtype=np.float64), np.ones( shape =(10), dtype=np.intc),np.ones( shape =(10), dtype=np.intc), np.ones( shape =(2,3), dtype=np.intc), None )

def test_expression_error_dim():
    assert_raises( ExpressionError, XTRAM, np.ones( shape =(2,3), dtype=np.intc), np.ones( shape =(2,10), dtype=np.float64 ), np.ones( shape =(10), dtype=np.intc ),np.ones( shape =(10), dtype=np.intc), np.ones( shape =(2,3), dtype=np.intc), np.ones( shape =(2), dtype=np.intc) )
    assert_raises( ExpressionError, XTRAM, np.ones( shape =(2,3,3), dtype=np.intc), np.ones( shape =(10), dtype=np.float64 ), np.ones( shape =(10), dtype=np.intc ),np.ones( shape =(10), dtype=np.intc), np.ones( shape =(2,3), dtype=np.intc), np.ones( shape =(2), dtype=np.intc) )
    assert_raises( ExpressionError, XTRAM, np.ones( shape =(2,3,3), dtype=np.intc), np.ones( shape =(2,10), dtype=np.float64 ), np.ones( shape =(10), dtype=np.intc ),np.ones( shape =(10), dtype=np.intc), np.ones( shape =(3), dtype=np.intc), np.ones( shape =(2), dtype=np.intc) )

def test_expression_error_markov():
    assert_raises( ExpressionError, XTRAM, np.ones( shape =(2,3,3), dtype=np.intc), np.ones( shape =(2,10), dtype=np.float64), np.ones( shape =(10), dtype=np.intc),np.ones( shape =(10), dtype=np.intc), np.ones( shape =(2,4), dtype=np.intc), np.ones( shape =(2), dtype=np.intc) )

def test_expression_error_therm():
    assert_raises( ExpressionError, XTRAM, np.ones( shape =(2,3,3), dtype=np.intc), np.ones( shape =(3,10), dtype=np.float64), np.ones( shape =(10), dtype=np.intc),np.ones( shape =(10), dtype=np.intc), np.ones( shape =(2,3), dtype=np.intc), np.ones( shape =(2), dtype=np.intc) )
    assert_raises( ExpressionError, XTRAM, np.ones( shape =(2,3,3), dtype=np.intc), np.ones( shape =(2,10), dtype=np.float64), np.ones( shape =(10), dtype=np.intc),np.ones( shape =(10), dtype=np.intc), np.ones( shape =(3,4), dtype=np.intc), np.ones( shape =(2), dtype=np.intc) )
    assert_raises( ExpressionError, XTRAM, np.ones( shape =(2,3,3), dtype=np.intc), np.ones( shape =(2,10), dtype=np.float64), np.ones( shape =(10), dtype=np.intc),np.ones( shape =(10), dtype=np.intc), np.ones( shape =(2,3), dtype=np.intc), np.ones( shape =(3), dtype=np.intc) )

'''def test_expression_error_int16():
    assert_raises( ExpressionError, XTRAM, np.ones( shape=(2,3,3), dtype=np.intc ), np.ones( shape=(2,3), dtype=np.int16 ) )

def test_expression_error_float32():
    assert_raises( ExpressionError, XTRAM, np.ones( shape=(2,3,3), dtype=np.intc ), np.ones( shape=(2,3), dtype=np.float32 ) )
    
    ( self, C_K_ij, u_I_x, T_x, M_x, N_K_i, N_K, target = 0, verbose = False ):
    
    assert_raises( ExpressionError, XTRAM, np.ones( shape =(2,3,3), dtype=np.intc), np.ones( shape =(3,10), dtype=np.float64), np.ones( shape =(10), dtype=np.intc),np.ones( shape =(10), dtype=np.intc) np.ones( shape =(3,3), dtype=np.intc), np.ones( shape =(3), dtype=np.intc) )
'''