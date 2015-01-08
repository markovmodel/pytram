.. _ref_api:

=================
API documentation
=================

.. toctree::
   :maxdepth: 2



The dTRAM estimator
===================

For running dTRAM, the following API functions are available on the package level, i.e., ::

   from pytram import dtram, dtram_me

dtram
-----

Run dTRAM using the TRAMData object as input:

.. autofunction:: pytram.api.dtram

dtram_me
--------

Run dTRAM using the mathematical expressions as input:

.. autofunction:: pytram.api.dtram_me

The xTRAM estimator
===================

For running xTRAM, the following API functions are available on the package level, i.e., ::

   from pytram import xtram, xtram_me

xtram
-----

Run xTRAM using the TRAMData object as input:

.. autofunction:: pytram.api.xtram

xtram_me
--------

Run xTRAM using the mathematical expressions as input:

.. autofunction:: pytram.api.xtram_me
