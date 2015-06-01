################################################################################
#
#   test_estimator.py - testing the basic estimator class
#
#   author: Christoph Wehmeyer <christoph.wehmeyer@fu-berlin.de>
#
################################################################################

from nose.tools import assert_raises, assert_equals
from pytram.estimator import Estimator, ExpressionError
import numpy as np

def test_expression_error_None():
    """test Estimator throws ExpressionError raise with None"""
    assert_raises(ExpressionError, Estimator, None)

def test_expression_error_int():
    """test Estimator throws ExpressionError raise with number"""
    assert_raises(ExpressionError, Estimator, 5)

def test_expression_error_list():
    """test Estimator throws ExpressionError raise with list"""
    assert_raises(ExpressionError, Estimator, [1, 2])

def test_expression_error_dim():
    """test Estimator throws ExpressionError raise with wrong dimension"""
    assert_raises(ExpressionError, Estimator, np.ones(shape=(2, 2), dtype=np.intc))

def test_expression_error_markov():
    """test Estimator throws ExpressionError raise with wrong Markov state count"""
    assert_raises(ExpressionError, Estimator, np.ones(shape=(2, 2, 3), dtype=np.intc))

def test_expression_error_float32():
    """test Estimator throws ExpressionError raise with wrong dtype"""
    assert_raises(ExpressionError, Estimator, np.ones(shape=(2, 3, 3), dtype=np.float32))

def test_expression_error_zeros():
    """test Estimator throws ExpressionError raise with zero counts"""
    assert_raises(ExpressionError, Estimator, np.zeros(shape=(2, 3, 3), dtype=np.intc))

def test_number_of_states():
    """test Estimator calculates state numbers"""
    estimator = Estimator(np.ones(shape=(2, 3, 3), dtype=np.intc))
    assert_equals(estimator.n_markov_states, 3)
    assert_equals(estimator.n_therm_states, 2)






