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
    I am the xTRAM estimator
    """
    def __init__( self, C_K_ij, u_I_x, T_x, M_x, N_K_i, N_K, target = 0, verbose = False ):

        r"""
        Initialize the XTRAM object
        
        Parameters
        ----------
        C_K_ij : 3-D numpy array
            Countmatrix for each thermodynamic state K
        u_I_x : 2-D numpy array
            Biasing tensor
        T_x : 1-D numpy array
            Thermodynamic state trajectory
        M_x : 1-D numpy array
            Markov state trajectories
        N_K_i : 2-D numpy array
            Number of markov samples in each thermodynamic state
        N_K : 1-D numpy array
            Numer of thermodynamic samples array
        target : Integer 
            target state for which pi_i should be computed
            default : 0
        verbose : Boolean
            Be loud and noisy
        """
        super( XTRAM, self ).__init__( C_K_ij )
        self._citation()
        self.verbose = verbose
        
        self.u_I_x = u_I_x
        self.T_x = T_x
        self.M_x = M_x
        self.N_K_i = N_K_i       
        self.N_K = N_K
        self.w_K = self._compute_w_K()
        self.f_K = self._compute_f_K()
        self.pi_K_i = self._compute_pi_K_i()
        self._ftol = 10e-15
        self._maxiter = 100000
        self.target = target
        self.pi_i = None
        
        
    def sc_iteration( self , ftol=10e-4, maxiter = 10, verbose = False):
        r"""Main iteration method
        Parameters
        ----------
            ftol : float
                tolerance of the free energy difference of each iteration update
                Default : 10e-10
            maxiter : int
                maximum number of iterations
                Default : 10
            verbose : boolean
                Be loud and noisy
                Default = False
        
        """
        finc = 0.0
        f_old = np.zeros(self.f_K.shape[0])
        self.b_i_IJ = np.zeros(shape=(self.n_markov_states, self.n_therm_states, self.n_therm_states))
       
        if verbose:
            print "# %8s %16s" % ( "[Step]", "[rel. Increment]" )
        for i in xrange( maxiter ):
            f_old[:]=self.f_K[:]
            b_i_IJ_equation( self.T_x, self.M_x, self.N_K, self.f_K, self.w_K, self.u_I_x, self.b_i_IJ )
            N_tilde = self._compute_sparse_N()
            C_i, C_j, C_ij, C_ji = self._compute_individual_N()
            x_row, c_column = self._initialise_X_and_N( N_tilde )
            ferr = iterate_x( N_tilde.shape[0],x_row.shape[0] , self._maxiter, self._ftol,C_i,C_j, C_ij, C_ji, x_row, c_column, x_row/x_row.sum() )
            #print 'ferr'+str( ferr )
            pi_curr =x_row/np.sum( x_row )
            self._update_pi_K_i( pi_curr )
            self._update_free_energies()
            finc = np.sum(np.abs( f_old-self.f_K) )
            #now we need to compute the results
            self.pi_i = self.pi_K_i[self.target]/self.pi_K_i[self.target].sum()
            if verbose:
                print "  %8d %16.8e" % ( i+1, finc )
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
                I_plus_one[I]+=self._metropolis(self.u_I_x[I,x],self.u_I_x[I+1,x])
            elif I==self.n_therm_states-1:
                I_minus_one[I]+=self._metropolis(self.u_I_x[I,x], self.u_I_x[I-1,x])
            else:
                I_plus_one[I]+=self._metropolis(self.u_I_x[I,x],self.u_I_x[I+1,x])
                I_minus_one[I]+=self._metropolis(self.u_I_x[I,x], self.u_I_x[I-1,x])
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
    # prints the needed citation                                       #
    #                                                                  #
    ####################################################################

    def _citation( self ):
        r"""Prints citation string"""
        citation_string = (
            "If you use this method for your data analysis please do not forget to cite:\n"
            "xTRAM: Estimating Equilibrium Expectations from Time-Correlated Simulation\n" 
            "Data at Multiple Thermodynamic States \n"
            "Antonia S. J. S. Mey, Hao Wu and Frank Noe, \n"
            "Phys. Rev. X 4, 041018")
        
        print citation_string

    ####################################################################
    #                                                                  #
    # Getters and setters and checks                                   #
    #                                                                  #
    ####################################################################

    @property
    def u_I_x( self ):
        return self._u_I_x
        
    @u_I_x.setter
    def u_I_x( self, u_I_x ):
        self._u_I_x = None
        if self._check_u_I_x( u_I_x ):
            self._u_I_x = u_I_x
    
    def _check_u_I_x( self, u_I_x ):
        if u_I_x is None:
            raise ExpressionError( "u_I_x", "is None" )
        if not isinstance( u_I_x, (np.ndarray,) ):
            raise ExpressionError( "u_I_x", "invalid type (%s)" % str( type( u_I_x ) ) )
        if 2 != u_I_x.ndim:
            raise ExpressionError( "u_I_x", "invalid number of dimensions (%d)" % u_I_x.ndim )
        if u_I_x.shape[0] != self.n_therm_states:
            raise ExpressionError( "u_I_x", "unmatching number of thermodynamic states (%d,%d)" % (u_I_x.shape[0], self.n_therm_states) )
        if np.float64 != u_I_x.dtype:
            raise ExpressionError( "u_I_x", "invalid dtype (%s)" % str( u_I_x.dtype ) )
        return True
    
    @property
    def M_x( self ):
        return self._M_x

    @M_x.setter
    def M_x( self, M_x ):
        self._M_x = None
        if self._check_M_x( M_x ):
            if self.verbose:
                print "M_x check pass"
            self._M_x = M_x

    def _check_M_x( self, M_x ):
        if M_x is None:
            raise ExpressionError( "M_x", "is None" )
        if not isinstance( M_x, (np.ndarray,) ):
            raise ExpressionError( "M_x", "invalid type (%s)" % str( type( M_x ) ) )
        if 1 != M_x.ndim:
            raise ExpressionError( "M_x", "invalid number of dimensions (%d)" % M_x.ndim )
        if M_x.shape[0] != self.u_I_x.shape[1]:
            raise ExpressionError( "M_x", "unmatching number thermodynamic samples (%d,%d)" % (M_x.shape[0], self.u_I_x.shape[1]) )
        if np.intc != M_x.dtype:
            raise ExpressionError( "M_x", "invalid dtype (%s)" % str( M_x.dtype ) )
        return True
        
    @property
    def T_x( self ):
        return self._T_x

    @T_x.setter
    def T_x( self, T_x ):
        self._T_x = None
        if self._check_T_x( T_x ):
            if self.verbose:
                print "T_x check pass"
            self._T_x = T_x

    def _check_T_x( self, T_x ):
        if T_x is None:
            raise ExpressionError( "T_x", "is None" )
        if not isinstance( T_x, (np.ndarray,) ):
            raise ExpressionError( "T_x", "invalid type (%s)" % str( type( T_x ) ) )
        if 1 != T_x.ndim:
            raise ExpressionError( "T_x", "invalid number of dimensions (%d)" % T_x.ndim )
        if T_x.shape[0] != self.u_I_x.shape[1]:
            raise ExpressionError( "T_x", "unmatching number thermodynamic samples (%d,%d)" % ( T_x.shape[0], self.u_I_x.shape[1] ) )
        if np.intc != T_x.dtype:
            raise ExpressionError( "T_x", "invalid dtype (%s)" % str( T_x.dtype ) )
        return True
        
    @property
    def N_K_i( self ):
        return self._N_K_i
        
    @N_K_i.setter
    def N_K_i( self, N_K_i ):
        self._N_K_i = None
        if self._check_N_K_i( N_K_i ):
            self._N_K_i = N_K_i
    
    def _check_N_K_i( self, N_K_i ):
        if N_K_i is None:
            raise ExpressionError( "N_K_i", "is None" )
        if not isinstance( N_K_i, (np.ndarray,) ):
            raise ExpressionError( "N_K_i", "invalid type (%s)" % str( type( N_K_i ) ) )
        if 2 != N_K_i.ndim:
            raise ExpressionError( "N_K_i", "invalid number of dimensions (%d)" % N_K_i.ndim )
        if N_K_i.shape[0] != self.n_therm_states:
            raise ExpressionError( "N_K_i", "unmatching number of thermodynamic states (%d,%d)" % ( N_K_i.shape[0], self.n_therm_states ) )
        if N_K_i.shape[1] != self.n_markov_states:
            raise ExpressionError( "N_K_i", "unmatching number of Markov states (%d,%d)" % ( N_K_i.shape[1], self.n_markov_states ) )
        return True
        
    @property
    def N_K( self ):
        return self._N_K
        
    @N_K.setter
    def N_K( self, N_K ):
        self._N_K = None
        if self._check_N_K( N_K ):
            self._N_K = N_K
    
    def _check_N_K( self, N_K ):
        if N_K is None:
            raise ExpressionError( "N_K", "is None" )
        if not isinstance( N_K, (np.ndarray,) ):
            raise ExpressionError( "N_K", "invalid type (%s)" % str( type( N_K ) ) )
        if 1 != N_K.ndim:
            raise ExpressionError( "N_K", "invalid number of dimensions (%d)" % N_K.ndim )
        if N_K.shape[0] != self.n_therm_states:
            raise ExpressionError( "N_K", "unmatching number of thermodynamic states (%d,%d)" % ( N_K.shape[0], self.n_therm_states) )
        return True
