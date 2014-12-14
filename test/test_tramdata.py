################################################################################
#
#   test_tramdata.py - testing the TRAMData class
#
#   author: Christoph Wehmeyer <christoph.wehmeyer@fu-berlin.de>
#
################################################################################

from nose.tools import assert_true
from pytram import TRAMData
import numpy as np

def test_Mx_Tx_nMTstates():
    seq = np.array( range( 20 ), dtype=np.intc )
    trajs = [ { 'm': seq[:10], 't': seq[:10] }, { 'm': seq[10:], 't': seq[10:] } ]
    td = TRAMData( trajs )
    assert_true( np.all( ( seq - td.T_x ) == 0  ) )
    assert_true( np.all( ( seq - td.M_x ) == 0  ) )
    assert_true( td.n_markov_states == 20 )
    assert_true( td.n_therm_states == 20 )

def test_N_K_i():
    trajs = []
    for K in xrange( 3 ):
        trajs.append( {
                'm': np.array( range( 4 ), dtype=np.intc ),
                't': np.ones( shape=(4,), dtype=np.intc )*K
            } )
    td = TRAMData( trajs )
    assert_true( np.all( td.N_K_i == 1 ) )
    assert_true( np.all( td.N_K == 4 ) )
    assert_true( td.N_K.sum() == 12 )
