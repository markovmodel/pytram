################################################################################
#
#   test_exceptions.py - testing the basic estimator class
#
#   author: Christoph Wehmeyer <christoph.wehmeyer@fu-berlin.de>
#
################################################################################

from nose.tools import assert_true
from pytram.estimator import ExpressionError, NotConvergedWarning

def test_expression_error():
    """test ExpressionError's attributes ad __str__"""
    ee = ExpressionError("Expression", "MSG")
    assert_true(ee.expression == "Expression")
    assert_true(ee.msg == "MSG")
    assert_true(ee.__str__() == "[Expression] MSG")

def test_not_converged_warning():
    """test NotConvergedWarning's attributes ad __str__"""
    ncw = NotConvergedWarning("Estimator", 0.1)
    assert_true(ncw.estimator == "Estimator")
    assert_true(ncw.increment == 0.1)
    assert_true(ncw.__str__() == "[Estimator] only reached increment %.3e" % 0.1)
