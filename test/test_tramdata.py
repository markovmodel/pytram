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
    """test M_x and T_x state assignments"""
    seq = np.array(range(20), dtype=np.intc)
    trajs = [{'m': seq[:10], 't': seq[:10]}, {'m': seq[10:], 't': seq[10:]}]
    td = TRAMData(trajs)
    assert_true(np.all((seq - td.T_x) == 0))
    assert_true(np.all((seq - td.M_x) == 0))
    assert_true(td.n_markov_states == 20)
    assert_true(td.n_therm_states == 20)

def test_N_K_i():
    """test states counts"""
    trajs = []
    for K in xrange(3):
        trajs.append({
            'm': np.array(range(4), dtype=np.intc),
            't': np.ones(shape=(4,), dtype=np.intc) * K})
    td = TRAMData(trajs)
    assert_true(np.all(td.N_K_i == 1))
    assert_true(np.all(td.N_K == 4))
    assert_true(td.N_K.sum() == 12)

def test_get_C_K_ij_lag1():
    """test transition counts at lag time 1"""
    trajs = [{'m': np.array([0, 1, 2, 0], dtype=np.intc), 't': np.zeros(shape=(4,), dtype=np.intc)}]
    td = TRAMData(trajs)
    C_K_ij = np.array([[[0, 1, 0], [0, 0, 1], [1, 0, 0]]], dtype=np.intc)
    assert_true(np.all(C_K_ij == td.get_C_K_ij(1)))

def test_get_C_K_ij_lag2():
    """test transition counts at lag time 2"""
    trajs = [{'m': np.array([0, 1, 2, 0], dtype=np.intc), 't': np.zeros(shape=(4,), dtype=np.intc)}]
    td = TRAMData(trajs)
    C_K_ij = np.array([[[0, 0, 1], [1, 0, 0], [0, 0, 0]]], dtype=np.intc)
    assert_true(np.all(C_K_ij == td.get_C_K_ij(2)))

def test_get_C_K_ij_lag3():
    """test transition counts at lag time 3"""
    trajs = [{'m': np.array([0, 1, 2, 0], dtype=np.intc), 't': np.zeros(shape=(4,), dtype=np.intc)}]
    td = TRAMData(trajs)
    C_K_ij = np.array([[[1, 0, 0], [0, 0, 0], [0, 0, 0]]], dtype=np.intc)
    assert_true(np.all(C_K_ij == td.get_C_K_ij(3)))
