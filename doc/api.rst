.. _ref_api:

=================
API documentation
=================

.. toctree::
   :maxdepth: 2



The dTRAM estimator
===================

For running dTRAM, the following API functions are available on the package level, i.e., ::

   from pytram import dtram, dtram_from_matrix

dtram
-----

Run dTRAM using the TRAMData object as input:

.. autofunction:: pytram.api.dtram

dtram_from_matrix
-----------------

Run dTRAM using the mathematical expressions as input:

.. autofunction:: pytram.api.dtram_from_matrix

The xTRAM estimator
===================

For running xTRAM, the following API functions are available on the package level, i.e., ::

   from pytram import xtram, xtram_from_matrix

xtram
-----

Run xTRAM using the TRAMData object as input:

.. autofunction:: pytram.api.xtram

xtram_from_matrix
-----------------

Run xTRAM using the mathematical expressions as input:

.. autofunction:: pytram.api.xtram_from_matrix
