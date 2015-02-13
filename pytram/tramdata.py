r"""
.. moduleauthor:: Antonia Mey <antonia.mey@fu-berlin.de>, Christoph Wehmeyer <christoph.wehmeyer@fu-berlin.de>

"""
from .estimator import ExpressionError
import numpy as np



####################################################################################################
#
#   TRAMDATA CLASS FOR STORING SEQUENTIAL SIMULATION DATA AND CONVERSION TO TRAM INPUT EXPRESSIONS
#
####################################################################################################

class TRAMData( object ):
    r"""
    I am the TRAMData class
    
    Notes
    -----
    I convert/process the trajectory list of dictionaries into data types that are useful
    
    """
    def __init__( self, trajs, b_K_i=None, kT_K=None, kT_target=None, verbose=False ):
        r"""
        Parameters
        ----------
        trajs : list of dictionaries
            each dictionary contains the following entries:
            'm' markov sequence in a 1-D numpy array of integers
            't' thermodynamic sequence in a 1-D numpy array of integers
            'u' reduced bias energy sequences in a 2-D numpy array of floats
        b_K_i : 2D numpy array 
            contains discrete reduced bias energies
            Default = None
        """
        self.trajs = trajs
        self._n_therm_states = None
        self._n_markov_states = None
        self._N_K_i = None
        self._N_K = None
        self._M_x = None
        self._T_x = None
        self._u_I_x = None
        self.b_K_i = b_K_i
        self.kT_K = kT_K
        self.kT_target = kT_target
        self.verbose=verbose

    ############################################################################
    #
    #   n_markov_states / n_therm_states getters
    #
    ############################################################################

    @property
    def n_markov_states( self ):
        if self._n_markov_states is None:
            if self.verbose:
                print "# Counting Markov states"
            self._n_markov_states = 0
            for traj in self.trajs:
                max_state = np.max( traj['m'] )
                if max_state > self._n_markov_states:
                    self._n_markov_states = max_state
            self._n_markov_states += 1
            if self.verbose:
                print "# ... found %d Markov states" % self._n_markov_states
        return self._n_markov_states

    @property
    def n_therm_states( self ):
        if self._n_therm_states is None:
            if self.verbose:
                print "# Counting thermodynamic states"
            self._n_therm_states = 0
            for traj in self.trajs:
                max_state = np.max( traj['t'] )
                if max_state > self._n_therm_states:
                    self._n_therm_states = max_state
            self._n_therm_states += 1
            if self.verbose:
                print "# ... found %d thermodynamic states" % self._n_therm_states
        return self._n_therm_states

    ############################################################################
    #
    #   N_K_i / N_K getters
    #
    ############################################################################

    @property
    def N_K_i( self ):
        if self._N_K_i is None:
            if self.verbose:
                print "# Counting visited Markov states"
            self._N_K_i = np.zeros( shape=(self.n_therm_states,self.n_markov_states), dtype=np.intc )
            for traj in self.trajs:
                for K in xrange( self.n_therm_states ):
                    inc_K = ( traj['t'] == K )
                    for i in xrange( self.n_markov_states ):
                        inc_i = ( traj['m'][inc_K] == i )
                        self._N_K_i[K,i] += inc_i.sum()
            if self.verbose:
                print "# ... done"
        return self._N_K_i

    @property
    def N_K( self ):
        if self._N_K is None:
            self._N_K = self.N_K_i.sum( axis=1 )
        return self._N_K.astype(np.intc)

    ############################################################################
    #
    #   M_x / T_x getters
    #
    ############################################################################

    @property
    def M_x( self ):
        if self._M_x is None:
            if self.verbose:
                print "# Copying Markov state sequences"
            self._M_x = np.zeros( shape=(self.N_K.sum(),), dtype=np.intc )
            a = 0
            for traj in self.trajs:
                b = a + traj['m'].shape[0]
                self._M_x[a:b] = traj['m'][:]
                a = b
            if self.verbose:
                print "# ... done"
        return self._M_x

    @property
    def T_x( self ):
        if self._T_x is None:
            if self.verbose:
                print "# Copying thermodynamic state sequences"
            self._T_x = np.zeros( shape=(self.N_K.sum(),), dtype=np.intc )
            a = 0
            for traj in self.trajs:
                b = a + traj['t'].shape[0]
                self._T_x[a:b] = traj['t'][:]
                a = b
            if self.verbose:
                print "# ... done"
        return self._T_x

    ############################################################################
    #
    #   C_K_ij getter method
    #
    ############################################################################

    def get_C_K_ij( self, lag, sliding_window=True ):
        C_K_ij = np.zeros(
                shape=(self.n_therm_states,self.n_markov_states,self.n_markov_states),
                dtype=np.intc
            )
        for traj in self.trajs:
            t = 0
            while t < traj['m'].shape[0]-lag:
                K = traj['t'][t]
                if np.all( traj['t'][t:t+lag+1] == K ):
                    C_K_ij[ K , traj['m'][t] , traj['m'][t+lag] ] += 1
                if sliding_window:
                    t += 1
                else:
                    t += lag
        return C_K_ij

    ############################################################################
    #
    #   u_I_x getter and helper functions
    #
    ############################################################################

    @property
    def u_I_x( self ):
        if self._u_I_x is None:
            if self.verbose:
                print "# Copying bias energy sequences"
            if ( None != self.kT_target ) and ( None != self.kT_K ):
                self.gen_u_I_x_from_kT_K()
            else:
                self.gen_u_I_x()
            if self.verbose:
                print "# ... done"
        return self._u_I_x

    def gen_u_I_x( self ):
        self._u_I_x = np.zeros( shape=(self.n_therm_states,self.N_K.sum()), dtype=np.float64 )
        a = 0
        for traj in self.trajs:
            b = a + traj['u'].shape[0]
            self._u_I_x[:,a:b] = traj['u'][:,:].transpose().copy()
            a = b

    def gen_u_I_x_from_kT_K( self ):
        u_x = np.zeros( shape=(self.N_K.sum(),), dtype=np.float64 )
        a = 0
        for traj in self.trajs:
            b = a + traj['u'].shape[0]
            u_x[a:b] = traj['u'][:,0]
            a = b
        for K in xrange( self.kT_K.shape[0] ):
            u_x[( self.T_x == K )] *= self.kT_K[K]
        self._u_I_x = np.zeros( shape=(self.kT_K.shape[0],self.N_K.sum()), dtype=np.float64 )
        for I in xrange( self.kT_K.shape[0] ):
            self._u_I_x[I,:] = ( 1.0/self.kT_K[I] - 1.0/self.kT_K[self.kT_target] ) * u_x[:]



















