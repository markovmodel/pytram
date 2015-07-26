################################################################################
#
#   test_three_state_model.py - testing pytram with a simple three state model
#
#   author: Christoph Wehmeyer <christoph.wehmeyer@fu-berlin.de>
#
################################################################################

from nose.tools import assert_raises, assert_true, assert_equal
from pytram import TRAMData, dtram, xtram
import numpy as np
from numpy.testing import assert_allclose

def tower_sample(distribution):
    """draws random integers from the given distribution"""
    cdf = np.cumsum(distribution)
    rnd = np.random.rand() * cdf[-1]
    ind = (cdf > rnd)
    idx = np.where(ind == True)
    return np.min(idx)

def evolve_chain(x, P, length):
    """generates a discrete Markov chain"""
    chain = np.zeros(length, dtype=np.intc)
    chain[0] = x
    for i in xrange(1, length):
        chain[i] = tower_sample(P[chain[i-1]])
    return chain

def assign_bias(dtraj, b_K_i):
    """assigns bias energies to discrete trajectories"""
    b = np.zeros(shape=(dtraj.shape[0], b_K_i.shape[0]), dtype=np.float64)
    for i in xrange(b_K_i.shape[1]):
        b[(dtraj == i), :] = (b_K_i[:, i])[np.newaxis, :]
    return b

def generate_data(P, b_K_i):
    """generates pyram compatible input data"""
    dtraj_0 = [evolve_chain(1, P[0, :, :], 100) for i in xrange(100)]
    dtraj_1 = [evolve_chain(1, P[1, :, :], 100) for i in xrange(100)]
    inp = [{
        'm': d,
        't': np.zeros(shape=d.shape, dtype=np.intc),
        'b': assign_bias(d, b_K_i)} for d in dtraj_0]
    inp += [{
        'm': d,
        't': np.ones(shape=d.shape, dtype=np.intc),
        'b': assign_bias(d, b_K_i)} for d in dtraj_1]
    return inp

class TestThreeStateModel(object):
    @classmethod
    def setup_class(cls):
        cls.energy = np.array([1.0, 2.0, 0.0], dtype=np.float64)
        cls.b_K_i = np.array([[0.0, 0.0, 0.0], 2.0 - cls.energy], dtype=np.float64)
        cls.pi_i = np.exp(-cls.energy) / np.exp(-cls.energy).sum()
        cls.f_i = -np.log(cls.pi_i)
        cls.F_K = 1.0 / (np.exp(-cls.b_K_i) * cls.pi_i[np.newaxis, :]).sum(axis=1)
        cls.f_K = np.log(cls.F_K)
        cls.pi_K_i = cls.F_K[:, np.newaxis] * np.exp(-cls.b_K_i) * cls.pi_i[np.newaxis, :]
        cls.f_K_i = -np.log(cls.pi_K_i)
        metropolis = cls.energy[np.newaxis, :] - cls.energy[:, np.newaxis]
        metropolis[(metropolis < 0.0)] = 0.0
        selection = np.array([[0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]], dtype=np.float64)
        metr_hast = selection * np.exp(-metropolis)
        for i in xrange(metr_hast.shape[0]):
            metr_hast[i, i] = 0.0
            metr_hast[i, i] = 1.0 - metr_hast[i, :].sum()
        cls.tmat = np.array([metr_hast, selection])
        cls.inp = generate_data(cls.tmat, cls.b_K_i)
    @classmethod
    def teardown_class(cls):
        pass
    def setup(self):
        pass
    def teardown(self):
        pass
    def test_dtram_api(self):
        """testing the dTRAM API"""
        tramdata = TRAMData(self.inp, b_K_i=self.b_K_i, verbose=True)
        dtram_obj = dtram(tramdata, lag=1, maxiter=1, ftol=1.0E-14, verbose=True)
        dtram_obj = dtram(tramdata, lag=1, maxiter=100000, ftol=1.0E-14, verbose=True)
        maxerr = 1.0E-1
        assert_allclose(dtram_obj.f_K, self.f_K, atol=maxerr)
        assert_allclose(dtram_obj.f_i, self.f_i, atol=maxerr)
        assert_allclose(dtram_obj.pi_i, self.pi_i, atol=maxerr)
        assert_allclose(dtram_obj.f_K_i, self.f_K_i, atol=maxerr)
        assert_allclose(dtram_obj.pi_K_i, self.pi_K_i, atol=maxerr)
        assert_allclose(dtram_obj.estimate_transition_matrices(), self.tmat, atol=maxerr)
    def test_xtram_api(self):
        """testing the xTRAM API"""
        tramdata = TRAMData(self.inp, verbose=True)
        xtram_obj = xtram(tramdata, lag=1, maxiter=1, ftol=1.0E-13, verbose=True)
        xtram_obj = xtram(tramdata, lag=1, maxiter=10000, ftol=1.0E-13, verbose=True)
        maxerr = 1.0E-1
        assert_allclose(xtram_obj.f_K, self.f_K, atol=maxerr)
        assert_allclose(xtram_obj.f_i, self.f_i, atol=maxerr)
        assert_allclose(xtram_obj.pi_i, self.pi_i, atol=maxerr)
        assert_allclose(xtram_obj.f_K_i, self.f_K_i, atol=maxerr)
        assert_allclose(xtram_obj.pi_K_i, self.pi_K_i, atol=maxerr)
