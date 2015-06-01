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
    """test ExpressionError raise"""
    assert_raises(ExpressionError, Estimator, None)

def test_expression_error_int():
    """test ExpressionError raise"""
    assert_raises(ExpressionError, Estimator, 5)

def test_expression_error_list():
    """test ExpressionError raise"""
    assert_raises(ExpressionError, Estimator, [1, 2])

def test_expression_error_dim():
    """test ExpressionError raise"""
    assert_raises(ExpressionError, Estimator, np.ones(shape=(2, 2), dtype=np.intc))

def test_expression_error_markov():
    """test ExpressionError raise"""
    assert_raises(ExpressionError, Estimator, np.ones(shape=(2, 2, 3), dtype=np.intc))

def test_expression_error_float32():
    """test ExpressionError raise"""
    assert_raises(ExpressionError, Estimator, np.ones(shape=(2, 3, 3), dtype=np.float32))

def test_expression_error_zeros():
    """test ExpressionError raise"""
    assert_raises(ExpressionError, Estimator, np.zeros(shape=(2, 3, 3), dtype=np.intc))

def test_number_of_states():
    """test calculation of state numbers"""
    estimator = Estimator(np.ones(shape=(2, 3, 3), dtype=np.intc))
    assert_equals(estimator.n_markov_states, 3)
    assert_equals(estimator.n_therm_states, 2)






