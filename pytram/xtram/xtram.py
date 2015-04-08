r"""

======================
xTRAM estimator module
======================

.. moduleauthor:: Antonia Mey <antonia.mey@fu-berlin.de>

"""

import numpy as np
from ..estimator import Estimator, NotConvergedWarning, ExpressionError
from .ext import b_i_IJ_equation, iterate_x




########################################################################
#                                                                      #
#   XTRAM ESTIMATOR CLASS                                              #
#                                                                      #
########################################################################


class XTRAM( Estimator ):
    r"""
    This is the XTRAM estimator class

    Parameters
    ----------
    C_K_ij : numpy.ndarray( shape=(T,M,M), dtype=numpy.intc )
        transition counts between the M discrete Markov states for each of the T thermodynamic ensembles
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
    def __init__( self, C_K_ij, b_K_x, T_x, M_x, N_K_i, target=0 ):

        super( XTRAM, self ).__init__( C_K_ij )
        if self._check_b_K_x( b_K_x ):
            self._b_K_x = b_K_x
        if self._check_M_x( M_x ):
            self._M_x = M_x
        if self._check_T_x( T_x ):
            self._T_x = T_x
        if self._check_N_K_i( N_K_i ):
            self._N_K_i = N_K_i.astype( np.intc )
        self._N_K = np.sum( N_K_i, axis=1 ).astype( np.intc )
        self.w_K = self._compute_w_K()
        self._f_K = self._compute_f_K()
        self._pi_K_i = self._compute_pi_K_i()
        self.target = target
        # citation information
        self.citation = [
                "xTRAM: Estimating Equilibrium Expectations from Time-Correlated Simulation",
                "Data at Multiple Thermodynamic States;",
                "Antonia S.J.S. Mey, Hao Wu, and Frank Noe",
                "Phys. Rev. X 4, 041018 (2014)"
            ]

    ############################################################################
    #
    #   override getters for stationary properties
    #
    ############################################################################

    @property
    def pi_i( self ):
        return self.pi_K_i[self.target] / self.pi_K_i[self.target].sum()

    @property
    def pi_K_i( self ):
        return self._pi_K_i

    @property
    def f_K( self ):
        return self._f_K

    ############################################################################
    #
    #   self-consistent-iteration
    #
    ############################################################################

    def sc_iteration( self , maxiter=100, ftol=1.0E-5, verbose = False):
        r"""
        sc_iteration function

        Parameters
        ----------
        maxiter : int (default=100)
            maximum number of self-consistent-iteration steps
        ftol : float (default=1.0E-5)
            convergence criterion based on the max relative change in an self-consistent-iteration step
        verbose : boolean (default=False)
            Be loud and noisy
        """
        finc = 0.0
        f_old = np.zeros( self.f_K.shape[0] )
        self.b_i_IJ = np.zeros( shape=(self.n_markov_states, self.n_therm_states, self.n_therm_states) )
        if verbose:
            print "# %25s %25s" % ( "[Step]", "[rel. Increment]" )
        for i in xrange( maxiter ):
            f_old[:]=self.f_K[:]
            b_i_IJ_equation( self.T_x, self.M_x, self.N_K, self.f_K, self.w_K, self.b_K_x, self.b_i_IJ )
            N_tilde = self._compute_sparse_N()
            C_i, C_j, C_ij, C_ji = self._compute_individual_N()
            x_row, c_column = self._initialise_X_and_N( N_tilde )
            ferr = iterate_x( N_tilde.shape[0], x_row.shape[0] , self._maxiter, self._ftol, C_i, C_j, C_ij, C_ji, x_row, c_column, x_row/x_row.sum() )
            #print 'ferr'+str( ferr )
            pi_curr = x_row / np.sum( x_row )
            self._update_pi_K_i( pi_curr )
            self._update_free_energies()
            finc = np.sum( np.abs( f_old - self.f_K ) )
            if verbose:
                print " %25d %25.12e" % ( i+1, finc )
            if finc < ftol:
                    break
        if finc > ftol:
                raise NotConvergedWarning( "XTRAM", finc )

    def _initialise_X_and_N( self, N_tilde ):
        r"""
            sets default values for x_i and N_i
        """
        X_row = np.zeros(np.max(N_tilde[:,0])+1)
        N_column = np.zeros(np.max(N_tilde[:,0])+1)
        for i in xrange(len(N_tilde)):
            entry = N_tilde[i]
            if entry[0]==entry[1]:
                X_row[entry[0]]+=(entry[2]+entry[3])*0.5
                N_column[entry[0]]+=entry[2]
            else:
                N_column[entry[0].astype(int)]+=entry[2] #Check that this is the right summation!
                N_column[entry[1].astype(int)]+=entry[3]
                X_row[entry[0]]+=(entry[2]+entry[3])*0.5
                X_row[entry[1]]+=(entry[2]+entry[3])*0.5
        return (X_row, N_column)

    def _update_pi_K_i( self, pi_curr ):
        r"""
        copies the current iteration pi_curr into the pi_K_i variable and normalises it as required
        """
        for K in xrange(self.n_therm_states):
            initial = K*self.n_markov_states
            final =K*self.n_markov_states+self.n_markov_states
            self.pi_K_i[K][:] = pi_curr[initial:final]/np.sum(pi_curr[:])

    def _update_free_energies( self ):
        r"""
        computes the free energies based on the current pi_K_i
        """
        for K in xrange( self.f_K.shape[0] ):
            self.f_K[K] = self.f_K[K]- np.log((np.sum(self.N_K).astype(float)/self.N_K[K])*(np.sum(self.pi_K_i[K,:])))

    ####################################################################
    #                                                                  #
    # Computes the extended count matrix                               #
    #                                                                  #
    ####################################################################

    def _compute_individual_N( self, factor=1.0 ):
        C_i = []
        C_j = []
        C_ij = []
        C_ji = []
        for I in xrange(self.n_therm_states):
            for i in xrange(self.n_markov_states):
                for j in xrange(i,self.n_markov_states):
                    s1=i+I*self.n_markov_states
                    s2=j+I*self.n_markov_states
                    if i==j:
                        n_ij = (self.C_K_ij[I,i,j]*factor+self.b_i_IJ[i,I,I])
                        n_ji = (self.C_K_ij[I,i,j]*factor+self.b_i_IJ[i,I,I])
                        C_i.append(s1)
                        C_j.append(s2)
                        C_ij.append(n_ij)
                        C_ji.append(n_ji)
                    else:
                        n_ij = self.C_K_ij[I,i,j]*factor
                        n_ji = self.C_K_ij[I,j,i]*factor
                        if n_ij or n_ji !=0: 
                           C_i.append(s1)
                           C_j.append(s2)
                           C_ij.append(n_ij)
                           C_ji.append(n_ji)
        for I in xrange(self.n_therm_states):
            for J in xrange(I,self.n_therm_states):
                for i in xrange(self.n_markov_states):
                    s1=self.n_markov_states*I+i
                    s2=self.n_markov_states*J+i
                    if I!=J:
                        n_ij = self.b_i_IJ[i,I,J]
                        n_ji = self.b_i_IJ[i,J,I]
                        C_i.append(s1)
                        C_j.append(s2)
                        C_ij.append(n_ij)
                        C_ji.append(n_ji)
        return (np.array(C_i).astype(np.intc), np.array(C_j).astype(dtype=np.intc), np.array(C_ij), np.array(C_ji))

    def _compute_sparse_N( self , factor=1.0):
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

        N_tilde=[]
        for I in xrange(self.n_therm_states):
            for i in xrange(self.n_markov_states):
                for j in xrange(i,self.n_markov_states):
                    s1=i+I*self.n_markov_states
                    s2=j+I*self.n_markov_states
                    if i==j:
                        n_ij = (self.C_K_ij[I,i,j]*factor+self.b_i_IJ[i,I,I])
                        n_ji = (self.C_K_ij[I,i,j]*factor+self.b_i_IJ[i,I,I])
                        entry=np.zeros(4)
                        entry[0] = s1
                        entry[1] = s2
                        entry[2] = n_ij
                        entry[3] = n_ji
                        N_tilde.append(entry)
                    else:
                        n_ij = self.C_K_ij[I,i,j]*factor
                        n_ji = self.C_K_ij[I,j,i]*factor
                        if n_ij or n_ji !=0: 
                            entry=np.zeros(4)
                            entry[0] = s1
                            entry[1] = s2
                            entry[2] = n_ij
                            entry[3] = n_ji
                            N_tilde.append(entry)
        for I in xrange(self.n_therm_states):
            for J in xrange(I,self.n_therm_states):
                for i in xrange(self.n_markov_states):
                    s1=self.n_markov_states*I+i
                    s2=self.n_markov_states*J+i
                    if I!=J:
                        n_ij = self.b_i_IJ[i,I,J]
                        n_ji = self.b_i_IJ[i,J,I]
                        entry=np.zeros(4)
                        entry[0] = s1
                        entry[1] = s2
                        entry[2] = n_ij
                        entry[3] = n_ji
                        N_tilde.append(entry)
        return np.array(N_tilde)

    ####################################################################
    #                                                                  #
    # Computes the initial guess of free energies vie bar ratios       #
    #                                                                  #
    ####################################################################

    def _compute_f_K( self ):
        _f_K = np.ones(self.n_therm_states)
        bar_ratio = self._bar_ratio()
        for I in xrange (1,self.n_therm_states):
            _f_K[I] = _f_K[I-1] - np.log(bar_ratio[I-1])
        return _f_K

    ####################################################################
    #                                                                  #
    # Computes BAR ratios                                              #
    #                                                                  #
    ####################################################################

    def _bar_ratio( self ):
        bar_ratio = np.zeros(self.n_therm_states-1)
        I_plus_one = np.zeros(self.n_therm_states)
        I_minus_one = np.zeros(self.n_therm_states)
        for x in xrange(self.T_x.shape[0]):
            I = self.T_x[x]
            if I==0:
                I_plus_one[I]+=self._metropolis(self.b_K_x[I,x],self.b_K_x[I+1,x])
            elif I==self.n_therm_states-1:
                I_minus_one[I]+=self._metropolis(self.b_K_x[I,x], self.b_K_x[I-1,x])
            else:
                I_plus_one[I]+=self._metropolis(self.b_K_x[I,x],self.b_K_x[I+1,x])
                I_minus_one[I]+=self._metropolis(self.b_K_x[I,x], self.b_K_x[I-1,x])
        for I in xrange(bar_ratio.shape[0]):
            bar_ratio[I]=(I_plus_one[I]/I_minus_one[I+1])*(self.N_K[I+1].astype('float')/self.N_K[I].astype('float'))
        return bar_ratio

    ####################################################################
    #                                                                  #
    # metropolis function                                              #
    #                                                                  #
    ####################################################################

    def _metropolis( self, u_1, u_2 ):
        if (u_1-u_2)>0:
            return 1.0
        else:
            return np.exp(u_1-u_2)

    ####################################################################
    #                                                                  #
    # Initialises the stationary probabilities                         #
    #                                                                  #
    ####################################################################
    
    def _compute_pi_K_i( self ):
        _pi_K_i = np.ones(self.n_therm_states*self.n_markov_states).reshape(self.n_therm_states,self.n_markov_states)
        return _pi_K_i

    ####################################################################
    #                                                                  #
    # Computes the the weight at each thermoydnamic state              #
    #                                                                  #
    ####################################################################
        
    def _compute_w_K( self ):
        return self.N_K.astype(float)/np.sum(self.N_K) #weight array based on thermodynamics sample counts

    ####################################################################
    #                                                                  #
    # Getters and setters and checks                                   #
    #                                                                  #
    ####################################################################

    @property
    def b_K_x( self ):
        return self._b_K_x

    def _check_b_K_x( self, b_K_x ):
        if b_K_x is None:
            raise ExpressionError( "b_K_x", "is None" )
        if not isinstance( b_K_x, (np.ndarray,) ):
            raise ExpressionError( "b_K_x", "invalid type (%s)" % str( type( b_K_x ) ) )
        if 2 != b_K_x.ndim:
            raise ExpressionError( "b_K_x", "invalid number of dimensions (%d)" % b_K_x.ndim )
        if b_K_x.shape[0] != self.n_therm_states:
            raise ExpressionError( "b_K_x", "not matching number of thermodynamic states (%d,%d)" % (b_K_x.shape[0], self.n_therm_states) )
        if np.float64 != b_K_x.dtype:
            raise ExpressionError( "b_K_x", "invalid dtype (%s)" % str( b_K_x.dtype ) )
        return True

    @property
    def M_x( self ):
        return self._M_x

    def _check_M_x( self, M_x ):
        if M_x is None:
            raise ExpressionError( "M_x", "is None" )
        if not isinstance( M_x, (np.ndarray,) ):
            raise ExpressionError( "M_x", "invalid type (%s)" % str( type( M_x ) ) )
        if 1 != M_x.ndim:
            raise ExpressionError( "M_x", "invalid number of dimensions (%d)" % M_x.ndim )
        if M_x.shape[0] != self.b_K_x.shape[1]:
            raise ExpressionError( "M_x", "not matching number thermodynamic samples (%d,%d)" % (M_x.shape[0], self.b_K_x.shape[1]) )
        if np.intc != M_x.dtype:
            raise ExpressionError( "M_x", "invalid dtype (%s)" % str( M_x.dtype ) )
        return True

    @property
    def T_x( self ):
        return self._T_x

    def _check_T_x( self, T_x ):
        if T_x is None:
            raise ExpressionError( "T_x", "is None" )
        if not isinstance( T_x, (np.ndarray,) ):
            raise ExpressionError( "T_x", "invalid type (%s)" % str( type( T_x ) ) )
        if 1 != T_x.ndim:
            raise ExpressionError( "T_x", "invalid number of dimensions (%d)" % T_x.ndim )
        if T_x.shape[0] != self.b_K_x.shape[1]:
            raise ExpressionError( "T_x", "not matching number thermodynamic samples (%d,%d)" % ( T_x.shape[0], self.b_K_x.shape[1] ) )
        if np.intc != T_x.dtype:
            raise ExpressionError( "T_x", "invalid dtype (%s)" % str( T_x.dtype ) )
        return True

    @property
    def N_K_i( self ):
        return self._N_K_i

    def _check_N_K_i( self, N_K_i ):
        if N_K_i is None:
            raise ExpressionError( "N_K_i", "is None" )
        if not isinstance( N_K_i, (np.ndarray,) ):
            raise ExpressionError( "N_K_i", "invalid type (%s)" % str( type( N_K_i ) ) )
        if 2 != N_K_i.ndim:
            raise ExpressionError( "N_K_i", "invalid number of dimensions (%d)" % N_K_i.ndim )
        if N_K_i.shape[0] != self.n_therm_states:
            raise ExpressionError( "N_K_i", "not matching number of thermodynamic states (%d,%d)" % ( N_K_i.shape[0], self.n_therm_states ) )
        if N_K_i.shape[1] != self.n_markov_states:
            raise ExpressionError( "N_K_i", "not matching number of Markov states (%d,%d)" % ( N_K_i.shape[1], self.n_markov_states ) )
        return True

    @property
    def N_K( self ):
        return self._N_K
