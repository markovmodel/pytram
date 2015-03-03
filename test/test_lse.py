################################################################################
#
#   test_lse.py - testing the logsumexp implementation
#
#   author: Christoph Wehmeyer <christoph.wehmeyer@fu-berlin.de>
#
################################################################################

from nose.tools import assert_raises, assert_true
from pytram.lse import sort, logsumexp, logsumexp_pair
import numpy as np

def test_sort_urnd():
    data = np.random.rand( 100000 ).astype( np.float64 )
    sort( data )
    for i in xrange( 1, data.shape[0] ):
        assert_true( data[i-1] <= data[i] )

def test_sort_grnd():
    data = np.random.randn( 100000 ).astype( np.float64 )
    sort( data )
    for i in xrange( 1, data.shape[0] ):
        assert_true( data[i-1] <= data[i] )

def test_logsumexp_converged_geometric_series():
    data = np.arange( 10000 ).astype( np.float64 )
    assert_true( np.abs( logsumexp( -data ) - 0.45867514538708193 ) < 1.0E-15 )
    
def test_logsumexp_truncated_diverging_gemoetric_series():
   data = np.arange( 10000 ).astype( np.float64 )
   assert_true( np.abs( logsumexp( data ) - 9999.4586751453862 ) < 1.0E-15 )

def test_logsumexp_pair():
    assert_true( np.abs( logsumexp_pair( 0.0, 0.0 ) - np.log( 2.0 ) ) < 1.0E-15 )
    assert_true( np.abs( logsumexp_pair( 1.0, 1.0 ) - ( 1.0 + np.log( 2.0 ) ) ) < 1.0E-15 )
    assert_true( np.abs( logsumexp_pair( 10.0, 10.0 ) - ( 10.0 + np.log( 2.0 ) ) ) < 1.0E-15 )
    assert_true( np.abs( logsumexp_pair( 100.0, 100.0 ) - ( 100.0 + np.log( 2.0 ) ) ) < 1.0E-15 )
    assert_true( np.abs( logsumexp_pair( 1000.0, 1000.0 ) - ( 1000.0 + np.log( 2.0 ) ) ) < 1.0E-15 )
    assert_true( np.abs( logsumexp_pair( 10.0, 0.0 ) - 10.000045398899218 ) < 1.0E-15 )
    assert_true( np.abs( logsumexp_pair( 0.0, 10.0 ) - 10.000045398899218 ) < 1.0E-15 )
    assert_true( np.abs( logsumexp_pair( 100.0, 0.0 ) - 100.0 ) < 1.0E-15 )
    assert_true( np.abs( logsumexp_pair( 0.0, 100.0 ) - 100.0 ) < 1.0E-15 )
    assert_true( np.abs( logsumexp_pair( 1000.0, 0.0 ) - 1000.0 ) < 1.0E-15 )
    assert_true( np.abs( logsumexp_pair( 0.0, 1000.0 ) - 1000.0 ) < 1.0E-15 )
