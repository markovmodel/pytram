#!/usr/bin/env python

####################################################################################################
#                                                                                                  #
#   RUN SCRIPT FOR THE DTRAM METHOD WITHIN THE PYTRAM package                                      #
#                                                                                                  #
#   author: Christoph Wehmeyer <christoph.wehmeyer@fu-berlin.de>                                   #
#                                                                                                  #
####################################################################################################



####################################################################################################
#
#   IMPORTS
#
####################################################################################################

from pytram import Reader, TRAMData, DTRAM, ExpressionError, NotConvergedWarning
from argparse import ArgumentParser, FileType
from sys import exit
import numpy as np



####################################################################################################
#
#   MAIN PART
#
####################################################################################################

if '__main__' == __name__:

    ############################################################################
    #
    #   capture the command line arguments
    #
    ############################################################################
    parser = ArgumentParser()
    parser.add_argument(
            'files',
            help='pytram compatible files for evaluation',
            nargs='*',
            metavar='FILE'
        )
    parser.add_argument(
            "--b_K_i_file",
            help="specify a pytram compatible file with discretised bias energies",
            metavar="FILE"
        )
    parser.add_argument(
            "--lag",
            help="specify a lag time for evaluation",
            type=int,
            default=1,
            metavar='INT'
        )
    parser.add_argument(
            "--maxlength",
            help="limit the number of trajectory frames",
            type=int,
            default=None,
            metavar='INT'
        )
    parser.add_argument(
            "--skiprows",
            help="Number of initial frames skipped",
            type=int,
            default=0,
            metavar='INT'
        )
    parser.add_argument(
            "--maxiter",
            help="limit the number of fixed point iterations",
            type=int,
            default=1000,
            metavar='INT'
        )
    parser.add_argument(
            "--ftol",
            help="limit the requested convergence level",
            type=float,
            default=1.0E-10,
            metavar='FLOAT'
        )
    parser.add_argument(
            "--verbose",
            help="show the progress during the self-consistent-iteration",
            action='store_true'
        )
    args = parser.parse_args()



    ############################################################################
    #
    #   check mandatory command line arguments
    #
    ############################################################################
    if args.b_K_i_file is None:
        print "ERROR: you must set the --b_K_i_file option!"
        exit( 1 )
    if 1 > len( args.files ):
        print "ERROR: you must give at least one pytram compatible trajectory file!"
        exit( 1 )


    ############################################################################
    #
    #   write header
    #
    ############################################################################
    print "\n\n###################################### PYTRAM ######################################"
    print "#\n#                          Invoking the dTRAM estimator"
    print "#\n### PARAMETERS\n#"
    print "# %25s %24d" % ( "[--lag]", args.lag )
    print "# %25s %24d" % ( "[--maxiter]", args.maxiter )
    print "# %25s %24.5e" % ( "[--ftol]", args.ftol )



    ############################################################################
    #
    #   import the data
    #
    ############################################################################
    print "#\n################################## IMPORTING DATA ##################################\n#"
    reader = Reader(
            args.files,
            b_K_i_file=args.b_K_i_file,
            maxlength=args.maxlength,
            skiprows=args.skiprows,
            verbose=True
        )
    tramdata = TRAMData( reader.trajs, b_K_i=reader.b_K_i )
    try:
        dtram_obj = DTRAM( tramdata.get_C_K_ij( args.lag ), tramdata.b_K_i )
    except ExpressionError, e:
        print "#\n### ERROR\n#"
        print "# Your input was faulty!"
        print "# The < %s > object is malformed: %s" % ( e.expression, e.msg )
        print "#\n### ABORTING\n\n"
        exit( 1 )
    print "#\n### SYSTEM INFORMATION\n#"
    print "# %25s %24d" % ( "[markov states]", tramdata.n_markov_states )
    print "# %25s %24d" % ( "[thermodynamic states]", tramdata.n_therm_states )



    ############################################################################
    #
    #   run the self-consistent-iteration
    #
    ############################################################################
    print "#\n#################################### RUN DTRAM #####################################\n#"
    try:
        print "# Run self-consistent-iteration"
        dtram_obj.sc_iteration( maxiter=args.maxiter, ftol=args.ftol, verbose=args.verbose )
        print "# ... converged!"
    except NotConvergedWarning, e:
        print "#\n### WARNING\n#\n# dTRAM is not converged - use these results carefuly!"
        print "#\n### RECOMMENDATION\n#\n# Run dtram.py again and increase --maxiter"



    ############################################################################
    #
    #   print out the results
    #
    ############################################################################
    print "#\n##################################### RESULTS ######################################"
    print "#\n### UNBIASED STATIONARY VECTOR\n#"
    print "# %25s %25s" % ( "[markov state]", "[stationary probability]" )
    for i in xrange( dtram_obj.pi_i.shape[0] ):
        print " %25d %25.12e" % ( i, dtram_obj.pi_i[i] )
    print "#\n### UNBIASED FREE ENERGY\n#"
    print "# %25s %25s" % ( "[markov state]", "[reduced free energy]" )
    f_i = -np.log( dtram_obj.pi_i )
    for i in xrange( f_i.shape[0] ):
        print " %25d %25.12e" % ( i, f_i[i] )
    print "#\n### THERMODYNAMIC FREE ENERGY\n#"
    print "# %25s %25s" % ( "[thermodynamic state]", "[reduced free energy]" )
    for i in xrange( dtram_obj.f_K.shape[0] ):
        print " %25d %25.12e" % ( i, dtram_obj.f_K[i] )



    ############################################################################
    #
    #   say good bye
    #
    ############################################################################
    print "#\n#################################### THAT'S IT #####################################\n#"
    print "#\n#                     Thank you for using the pytram package!\n#\n#"
    print "### CITATION\n#"
    print "#    Statistically optimal analysis of state-discretized trajectory data"
    print "#    from multiple thermodynamic states;"
    print "#    Hao Wu, Antonia S.J.S. Mey, Edina Rosta, and Frank Noe"
    print "#    J. Chem. Phys. 141, 214106 (2014)"
    print "#\n####################################################################################\n\n"











