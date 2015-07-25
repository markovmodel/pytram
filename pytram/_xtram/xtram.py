r"""

======================
xTRAM estimator module
======================

.. moduleauthor:: Antonia Mey <antonia.mey@fu-berlin.de>

"""

import numpy as np
import warnings
from pytram.estimator import Estimator, NotConvergedWarning, ExpressionError
from .ext import b_i_IJ_equation, iterate_x




########################################################################
#                                                                      #
#   XTRAM ESTIMATOR CLASS                                              #
#                                                                      #
########################################################################


class XTRAM(Estimator):
    r"""
    This is the XTRAM estimator class

    Parameters
    ----------
    C_K_ij : numpy.ndarray( shape=(T,M,M), dtype=numpy.intc )
        transition counts between the M discrete Markov states for each of the T
        thermodynamic ensembles
    b_K_x : numpy.ndarray( shape=(T,nsamples), dtype=float )
        Biasing tensor
    T_x : numpy.ndarray( shape=(nsamples), dtype=float )
        Thermodynamic state trajectory
    M_x : numpy.ndarray( shape=(nsamples), dtype=numpy.intc )
        Markov state trajectories
    N_K_i : numpy.ndarray( shape=(T,M), dtype=float )
        Number of markov samples (M) in each thermodynamic state (T)
    target : int (default=0)
        target state for which pi_i should be computed
    """
    def __init__(self, C_K_ij, b_K_x, T_x, M_x, N_K_i, target=0):

        super(XTRAM, self).__init__(C_K_ij)
        if self._check_b_K_x(b_K_x):
            self._b_K_x = b_K_x
        if self._check_M_x(M_x):
            self._M_x = M_x
        if self._check_T_x(T_x):
            self._T_x = T_x
        if self._check_N_K_i(N_K_i):
            self._N_K_i = N_K_i.astype(np.intc)
        self._N_K = np.sum(N_K_i, axis=1).astype(np.intc)
        self._N = np.sum(self._N_K)
        self.n_samples = len(M_x)
        self.w_K = self._compute_w_K()
        self._f_K = self._init_f_K()
        self._pi_K_i = self._compute_pi_K_i()
        self.target = target
        # inner iteration
        self._maxiter_inner = 100000
        self._ftol_inner = 1.0e-15
        # citation information
        self.citation = [
            "xTRAM: Estimating Equilibrium Expectations from Time-Correlated Simulation",
            "Data at Multiple Thermodynamic States;",
            "Antonia S.J.S. Mey, Hao Wu, and Frank Noe",
            "Phys. Rev. X 4, 041018 (2014)"]

    ############################################################################
    #
    #   override getters for stationary properties
    #
    ############################################################################

    @property
    def pi_i(self):
        return self.pi_K_i[self.target]

    @property
    def pi_K_i(self):
        return self._pi_K_i / self._pi_K_i.sum(axis=1)[:, np.newaxis]

    @property
    def f_K(self):
        return self._f_K

    ############################################################################
    #
    #   self-consistent-iteration
    #
    ############################################################################

    def sc_iteration(self, maxiter=100, ftol=1.0E-5, verbose=False):
        r"""
        sc_iteration function

        Parameters
        ----------
        maxiter : int (default=100)
            maximum number of self-consistent-iteration steps
        ftol : float (default=1.0E-5)
            convergence criterion based on the max relative change in a
            self-consistent-iteration step
        verbose : boolean (default=False)
            Be loud and noisy
        """
        finc = 0.0
        f_old = np.zeros(self.f_K.shape[0])
        self.b_i_IJ = np.zeros((self.n_markov_states, self.n_therm_states, self.n_therm_states), dtype=np.float64)
        if verbose:
            print "# %25s %25s" % ("[Step]", "[rel. Increment]")
        for i in xrange(maxiter):
            f_old[:] = self.f_K[:]
            # compute thermodynamic count matrix
            self.b_i_IJ = self._count_matrices_thermo()
            # TODO: should be fixed in C code that is currently skipped
            # b_i_IJ_equation(
            #     self.T_x, self.M_x, self.N_K, self.f_K, self.w_K, self.b_K_x, self.b_i_IJ)
            N_tilde = self._compute_sparse_N()
            C_i, C_j, C_ij, C_ji = self._compute_individual_N()
            x_row, c_column = self._initialise_X_and_N(N_tilde)
            ferr = iterate_x(
                N_tilde.shape[0],
                x_row.shape[0],
                self._maxiter_inner,
                self._ftol_inner,
                C_i,
                C_j,
                C_ij,
                C_ji,
                x_row,
                c_column,
                x_row/x_row.sum())
            #print 'ferr'+str( ferr )
            pi_curr = x_row / np.sum(x_row)
            self._update_pi_K_i(pi_curr)
            self._update_free_energies()
            finc = np.sum(np.abs(f_old - self.f_K))
            if verbose:
                print " %25d %25.12e" % (i+1, finc)
            if finc < ftol:
                break
        if finc > ftol:
            warnings.warn("XTRAM only reached increment %.3e" % finc)

    def _initialise_X_and_N(self, N_tilde):
        r"""
            sets default values for x_i and N_i
        """
        X_row = np.zeros(int(np.max(N_tilde[:, 0])+1))
        N_column = np.zeros(int(np.max(N_tilde[:, 0])+1))
        for i in xrange(len(N_tilde)):
            entry = N_tilde[i]
            if entry[0] == entry[1]:
                X_row[int(entry[0])] += (entry[2] + entry[3]) * 0.5
                N_column[int(entry[0])] += entry[2]
            else:
                N_column[int(entry[0])] += entry[2] #Check that this is the right summation!
                N_column[int(entry[1])] += entry[3]
                X_row[int(entry[0])] += (entry[2] + entry[3]) * 0.5
                X_row[int(entry[1])] += (entry[2] + entry[3]) * 0.5
        return (X_row, N_column)

    def _update_pi_K_i(self, pi_curr):
        r"""
        copies the current iteration pi_curr into the pi_K_i variable and normalises it as required
        """
        for K in xrange(self.n_therm_states):
            initial = K * self.n_markov_states
            final = K * self.n_markov_states + self.n_markov_states
            self._pi_K_i[K][:] = pi_curr[initial:final] / np.sum(pi_curr[:])

    def _update_free_energies(self):
        r"""
        computes the free energies based on the current pi_K_i
        """
        for K in xrange(self.f_K.shape[0]):
            self._f_K[K] = self._f_K[K] - np.log(
                (np.sum(self.N_K).astype(float) / self.N_K[K]) * (np.sum(self._pi_K_i[K, :])))

    ####################################################################
    #                                                                  #
    # Computes the extended count matrix                               #
    #                                                                  #
    ####################################################################

    #def _count_matrices_conf(ttrajs, dtrajs, lag):
    #    import msmtools.estimation as msmest
    #    nthermo = msmest.number_of_states(ttrajs)
    #    nstates = msmest.number_of_states(dtrajs)
    #    ntrajs = len(ttrajs)
    #    Cs = np.zeros((nthermo, nstates, nstates), dtype=np.intc)
    #    for i in xrange(ntrajs):
    #        ttraj = ttrajs[i]
    #        dtraj = dtrajs[i]
    #        for t in xrange(len(ttraj)-lag):
    #            Cs[ttraj[t], dtraj[t], dtraj[t+lag]] += 1
    #    return Cs

    def _count_matrices_thermo(self, metropolis=True):
        """Computes count matrix between thermodynamic states
        """
        N_k = self._N_K.astype(float)
        N = np.sum(N_k)
        Bs = np.zeros((self.n_markov_states, self.n_therm_states, self.n_therm_states), dtype=np.float64)
        for I in range(self.n_therm_states):
            # get samples starting from I
            indI = np.where(self._T_x == I)[0]
            # look at all targets
            p_IJ = np.zeros((self.n_therm_states, len(indI)))
            for J in range(self.n_therm_states):
                if I != J:
                    if metropolis:
                        p_IJ[J] = np.minimum(1.0, (N_k[J]/N_k[I]) * np.exp(self._f_K[J] - self._b_K_x[J, indI]
                                                                           - self._f_K[I] + self._b_K_x[I, indI]))
                        p_IJ[J] *= N_k[I]/N
                    else:
                        raise NotImplementedError()
            p_IJ[I] = np.ones(len(indI)) - p_IJ.sum(axis=0)
            # accumulate counts by discrete state
            d_arr_i = self._M_x[indI]
            for i in range(self.n_markov_states):
                indi = np.where(d_arr_i == i)[0]
                Bs[i, I, :] = p_IJ[:, indi].sum(axis=1)
        return Bs

    def _compute_individual_N(self, factor=1.0):
        C_i = []
        C_j = []
        C_ij = []
        C_ji = []
        for I in xrange(self.n_therm_states):
            for i in xrange(self.n_markov_states):
                for j in xrange(i,self.n_markov_states):
                    s1 = i + I * self.n_markov_states
                    s2 = j + I * self.n_markov_states
                    if i == j:
                        n_ij = (self.C_K_ij[I, i, j] * factor + self.b_i_IJ[i, I, I])
                        n_ji = (self.C_K_ij[I, i, j] * factor + self.b_i_IJ[i, I, I])
                        C_i.append(s1)
                        C_j.append(s2)
                        C_ij.append(n_ij)
                        C_ji.append(n_ji)
                    else:
                        n_ij = self.C_K_ij[I, i, j] * factor
                        n_ji = self.C_K_ij[I, j, i] * factor
                        if n_ij or n_ji != 0: # TODO: CHECK THIS LINE
                           C_i.append(s1)
                           C_j.append(s2)
                           C_ij.append(n_ij)
                           C_ji.append(n_ji)
        for I in xrange(self.n_therm_states):
            for J in xrange(I, self.n_therm_states):
                for i in xrange(self.n_markov_states):
                    s1 = self.n_markov_states * I + i
                    s2 = self.n_markov_states * J + i
                    if I != J:
                        n_ij = self.b_i_IJ[i, I, J]
                        n_ji = self.b_i_IJ[i, J, I]
                        C_i.append(s1)
                        C_j.append(s2)
                        C_ij.append(n_ij)
                        C_ji.append(n_ji)
        return (
            np.array(C_i).astype(np.intc),
            np.array(C_j).astype(dtype=np.intc),
            np.array(C_ij),
            np.array(C_ji))

    def _compute_sparse_N(self , factor=1.0):
        r"""Computes a Nx4 array containing the count matrix in a sparse format
        
        Parameters
        ----------
        factor : float
            multiplication factor default of 1 is fine

        Returns
        -------
        N_tilde : numpy 2d-array
            N-4 numpy array containing the count matrix N-tilde
        """
        N_tilde = []
        for I in xrange(self.n_therm_states):
            for i in xrange(self.n_markov_states):
                for j in xrange(i, self.n_markov_states):
                    s1 = i + I * self.n_markov_states
                    s2 = j + I * self.n_markov_states
                    if i == j:
                        n_ij = (self.C_K_ij[I, i, j] * factor + self.b_i_IJ[i, I, I])
                        n_ji = (self.C_K_ij[I, i, j] * factor + self.b_i_IJ[i, I, I])
                        entry = np.zeros(4)
                        entry[0] = s1
                        entry[1] = s2
                        entry[2] = n_ij
                        entry[3] = n_ji
                        N_tilde.append(entry)
                    else:
                        n_ij = self.C_K_ij[I, i, j] * factor
                        n_ji = self.C_K_ij[I, j, i] * factor
                        if n_ij or n_ji != 0: # TODO: CHECK THIS LINE
                            entry = np.zeros(4)
                            entry[0] = s1
                            entry[1] = s2
                            entry[2] = n_ij
                            entry[3] = n_ji
                            N_tilde.append(entry)
        for I in xrange(self.n_therm_states):
            for J in xrange(I, self.n_therm_states):
                for i in xrange(self.n_markov_states):
                    s1 = self.n_markov_states * I + i
                    s2 = self.n_markov_states * J + i
                    if I != J:
                        n_ij = self.b_i_IJ[i, I, J]
                        n_ji = self.b_i_IJ[i, J, I]
                        entry = np.zeros(4)
                        entry[0] = s1
                        entry[1] = s2
                        entry[2] = n_ij
                        entry[3] = n_ji
                        N_tilde.append(entry)
        return np.array(N_tilde)


    def _init_f_K(self):
        """ Computes the initial guess of free energies via bar ratios
        """
        I_plus_one = np.zeros(self.n_therm_states)
        I_minus_one = np.zeros(self.n_therm_states)
        for x in xrange(self.n_samples):
            I = self._T_x[x]
            if I > 0:
                I_minus_one[I] += self._metropolis(self._b_K_x[I, x], self._b_K_x[I-1, x])
            if I < self.n_therm_states-1:
                I_plus_one[I] += self._metropolis(self._b_K_x[I, x], self._b_K_x[I+1, x])
        # compute BAR free energies
        f_K = np.zeros(self.n_therm_states)
        for I in xrange(1, self.n_therm_states):
            bar_ratio = (I_plus_one[I-1] / float(self._N_K[I-1])) / (I_minus_one[I] / float(self._N_K[I]))
            f_K[I] = f_K[I-1] - np.log(bar_ratio)
        return f_K

    def _metropolis(self, u_1, u_2):
        """ Metropolis function
        """
        if (u_1 - u_2) > 0:
            return 1.0
        else:
            return np.exp(u_1 - u_2)

    def _compute_pi_K_i(self):
        """Initializes the stationary probabilities
        """
        _pi_K_i = np.ones(
            self.n_therm_states * self.n_markov_states).reshape(
                self.n_therm_states, self.n_markov_states)
        return _pi_K_i

    def _compute_w_K(self):
        """Computes the the weight at each thermodynamic state """
        return self.N_K.astype(np.float64) / np.sum(self.N_K)
        #weight array based on thermodynamics sample counts

    ####################################################################
    #                                                                  #
    # Getters and setters and checks                                   #
    #                                                                  #
    ####################################################################

    @property
    def b_K_x(self):
        return self._b_K_x

    def _check_b_K_x(self, b_K_x):
        if b_K_x is None:
            raise ExpressionError("b_K_x", "is None")
        if not isinstance(b_K_x, (np.ndarray,)):
            raise ExpressionError("b_K_x", "invalid type (%s)" % str(type(b_K_x)))
        if 2 != b_K_x.ndim:
            raise ExpressionError("b_K_x", "invalid number of dimensions (%d)" % b_K_x.ndim)
        if b_K_x.shape[0] != self.n_therm_states:
            raise ExpressionError("b_K_x", "not matching number of thermodynamic states (%d,%d)" \
                % (b_K_x.shape[0], self.n_therm_states))
        if np.float64 != b_K_x.dtype:
            raise ExpressionError("b_K_x", "invalid dtype (%s)" % str(b_K_x.dtype))
        return True

    @property
    def M_x(self):
        return self._M_x

    def _check_M_x(self, M_x):
        if M_x is None:
            raise ExpressionError("M_x", "is None")
        if not isinstance(M_x, (np.ndarray,)):
            raise ExpressionError("M_x", "invalid type (%s)" % str(type(M_x)))
        if 1 != M_x.ndim:
            raise ExpressionError("M_x", "invalid number of dimensions (%d)" % M_x.ndim)
        if M_x.shape[0] != self.b_K_x.shape[1]:
            raise ExpressionError("M_x", "not matching number thermodynamic samples (%d,%d)" \
                % (M_x.shape[0], self.b_K_x.shape[1]))
        if np.intc != M_x.dtype:
            raise ExpressionError("M_x", "invalid dtype (%s)" % str(M_x.dtype))
        return True

    @property
    def T_x(self):
        return self._T_x

    def _check_T_x(self, T_x):
        if T_x is None:
            raise ExpressionError("T_x", "is None")
        if not isinstance(T_x, (np.ndarray,)):
            raise ExpressionError("T_x", "invalid type (%s)" % str(type(T_x)))
        if 1 != T_x.ndim:
            raise ExpressionError("T_x", "invalid number of dimensions (%d)" % T_x.ndim)
        if T_x.shape[0] != self.b_K_x.shape[1]:
            raise ExpressionError("T_x", "not matching number thermodynamic samples (%d,%d)" \
                % (T_x.shape[0], self.b_K_x.shape[1]))
        if np.intc != T_x.dtype:
            raise ExpressionError("T_x", "invalid dtype (%s)" % str(T_x.dtype))
        return True

    @property
    def N_K_i(self):
        return self._N_K_i

    def _check_N_K_i(self, N_K_i):
        if N_K_i is None:
            raise ExpressionError("N_K_i", "is None")
        if not isinstance(N_K_i, (np.ndarray,)):
            raise ExpressionError("N_K_i", "invalid type (%s)" % str(type(N_K_i)))
        if 2 != N_K_i.ndim:
            raise ExpressionError("N_K_i", "invalid number of dimensions (%d)" % N_K_i.ndim)
        if N_K_i.shape[0] != self.n_therm_states:
            raise ExpressionError("N_K_i", "not matching number of thermodynamic states (%d,%d)" \
                % (N_K_i.shape[0], self.n_therm_states))
        if N_K_i.shape[1] != self.n_markov_states:
            raise ExpressionError("N_K_i", "not matching number of Markov states (%d,%d)" \
                % (N_K_i.shape[1], self.n_markov_states))
        return True

    @property
    def N_K(self):
        return self._N_K
