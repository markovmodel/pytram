.. _ref_user:

==========
User guide
==========

.. toctree::
   :maxdepth: 2

Running dTRAM
=============

The dTRAM estimator requires count matrices for all thermodynamic states (`C_K_ij`) and discrete reduced bias energies (`b_K_i`). ::

   # import the dTRAM API function
   from pytram import dtram_me

   # create the dTRAM object
   dtram_obj = dtram_me( C_K_ij, b_K_i )

   # show unbiased stationary distribution
   print dtram_obj.pi_i

   # get transition matrix for thermodynamic state 0
   T_0 = dtram_obj.estimate_transtition_matrix( 0 )
   print T_0

   # show thermodynamic free energies
   print dtram_obj.f_K

The difficult part lies in preparing the count matrices. This, however, can be taken care of by pytram...



