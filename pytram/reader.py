r"""
.. moduleauthor:: Christoph Wehmeyer <christoph.wehmeyer@fu-berlin.de>, Antonia Mey <antonia.mey@fu-berlin.de>

"""
import numpy as np



####################################################################################################
#
#   READER CLASS FOR IMPORTING SEQUENTIAL SIMULATION DATA
#
####################################################################################################

class Reader( object ):
    r"""
    Parameters
    ----------
    files : array_like
        list of filenames of the to-be-imported trajectory files
    b_K_i_file : string (optional)
        name of the file with the discretised reduced bias energies (b_K_i)
    kT_file : string (optional)
        name of the file with kT values from a multi-temperature simulation
    skiprows : int (default=0)
        skip the leading lines
    maxlength : int (optional)
        limit the maximal number of samples to use
    verbose : boolean (default=False)
        show import progress

    Notes
    -----
    I import trajectories from a list of files
    """
    def __init__( self, files, b_K_i_file=None, kT_file=None, skiprows=0, maxlength=None, verbose=False ):
        self.files = files
        self.b_K_i_file = b_K_i_file
        self.kT_file = kT_file
        self.maxlength = maxlength
        self.skiprows = skiprows
        self.verbose = verbose
        self._trajs = None
        self._b_K_i = None
        self._kT_K= None

    ############################################################################
    #
    #   trajs getter
    #
    ############################################################################

    @property
    def trajs( self ):
        if self._trajs is None:
            self._trajs = []
            for f in self.files:
                if self.verbose:
                    print "# Reading file <%s>" % f
                try:
                    content = np.loadtxt( f, dtype=np.float64, skiprows=self.skiprows )
                except IOError, e:
                    print "# ... cannot read file <%s> (ignored)" % f
                    continue
                length = content.shape[0]
                if self.verbose:
                    print "# ... length=%d" % length
                if ( None is not self.maxlength ) and ( self.maxlength < length ):
                    length = self.maxlength
                    if self.verbose:
                        print "# ... truncating to length=%d" % self.maxlength
                m = content[:length,0].astype( np.intc, copy=True )
                t = content[:length,1].astype( np.intc, copy=True )
                b = None
                if content.shape[1] > 2:
                    b = np.copy( content[:length,2:] )
                if self.verbose:
                    if None == b:
                        print "# ... no energy data"
                    else:
                        print "# ... %d energy column(s)" % ( content.shape[1] - 2 )
                self._trajs.append( { 'm': m, 't': t, 'b': b } )
        return self._trajs

    ############################################################################
    #
    #   b_K_i getter
    #
    ############################################################################

    @property
    def b_K_i( self ):
        if self._b_K_i is None:
            if self.b_K_i_file is None:
                return None
            if self.verbose:
                print "# Reading b_K_i_file <%s>" % self.b_K_i_file
            try:
                self._b_K_i = np.loadtxt( self.b_K_i_file, dtype=np.float64 ).transpose().copy()
                if self.verbose:
                    print "# ... found %d markov and %d thermodynamic states" % ( self._b_K_i.shape[1], self._b_K_i.shape[0] )
            except IOError, e:
                print "# ... cannot read file <%s>" % self.b_K_i_file
        return self._b_K_i

    ############################################################################
    #
    #   kT getter
    #
    ############################################################################

    @property
    def kT_K( self ):
        if self._kT_K is None:
            if self.kT_file is None:
                return None
            if self.verbose:
                print "# Reading kT_file <%s>" % self.kT_file
            try:
                self._kT_K = np.loadtxt( self.kT_file, dtype=np.float64 )
                if self._kT_K.ndim > 1:
                    if self.verbose:
                        print "# ... found %d columns - restricting to first column" % self._kT_K.shape[1]
                    self._kT_K = self._kT_K[:,0].copy()
                if self.verbose:
                    print "# ... found %d thermodynamic states" % self._kT_K.shape[0]
            except IOError, e:
                print "# ... cannot read file <%s>" % self.kT_file
        return self._kT_K


