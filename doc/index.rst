.. pytram documentation master file, created by
   sphinx-quickstart on Sat Nov 29 23:53:07 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

==================
The pytram package
==================

This software offers the python interface to the TRAM methods family. TRAM (transition-based reweighting analysis method) is a collection of Markov state model (MSM) estimators for analysing multi-ensemble simulations. The pytram package is implemented mostly in `NumPy <http://www.numpy.org/>`_ and `SciPy <http://www.scipy.org>`_, with `Cython <http://www.cython.org/>`_-based C-extensions for numerically demanding tasks. For documentation of the API, please have a look at the :ref:`ref_api`. To install this software and additional dependencies refer to the :ref:`Installation Guide <ref_install>`. The :ref:`User Guide <ref_user>` is a good starting point to learn how to use TRAM. For support/bug reports/sugguestions/complains, please visit us at `GitHub <http://github.com/markovmodel/pytram/>`_.

Contents:

.. toctree::
   :maxdepth: 2

   install
   user
   api
   reader
   tramdata
