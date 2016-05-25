******
pytram
******

The pytram package is decprecated and no longer supported. We recommend to switch to `PyEMMA <https://github.com/markovmodel/PyEMMA>`_.

.. image:: http://badges.github.io/stability-badges/dist/deprecated.svg
   :target: http://github.com/badges/stability-badges
.. image:: https://travis-ci.org/markovmodel/pytram.svg?branch=devel
   :target: https://travis-ci.org/markovmodel/pytram
.. image:: https://coveralls.io/repos/markovmodel/pytram/badge.svg?branch=devel
   :target: https://coveralls.io/r/markovmodel/pytram?branch=devel
.. image:: https://badge.fury.io/py/pytram.svg
   :target: https://pypi.python.org/pypi/pytram

This python package implements the transition-based reweighting analyis method (TRAM) estimators.



Installation
============

Using conda::

   conda install -c https://conda.binstar.org/omnia pytram


Using pip from PyPI::

   # you might have to install these dependencies manually
   pip install cython
   pip install numpy

   # install pytram
   pip install pytram


Using pip from github (this will install the latest development version)::

   # you might have to install these dependencies manually
   pip install cython
   pip install numpy

   # install pytram - this might be slow
   pip install git+https://github.com/markovmodel/pytram.git@devel


Authors
=======

- Christoph Wehmeyer <christoph.wehmeyer@fu-berlin.de>
- Antonia Mey <antonia.mey@fu-berlin.de>
- Fabian Paul <fab@physik.tu-berlin.de>
- Hao Wu <hwu@zedat.fu-berlin.de>
- Frank Noé <frank.noe@fu-berlin.de>



References
==========

* **dTRAM**: *Statistically optimal analysis of state-discretized trajectory data from multiple thermodynamic states*, Hao Wu, Antonia S.J.S. Mey, Edina Rosta, and Frank Noé, **J. Chem. Phys.** 141, 214106 (2014). 

    Download: <http://scitation.aip.org/content/aip/journal/jcp/141/21/10.1063/1.4902240>

* **xTRAM**: *Estimating Equilibrium Expectations from Time-Correlated Simulation Data at Multiple Thermodynamic States*, Antonia S.J.S. Mey, Hao Wu, and Frank Noé, **Phys. Rev. X** 4, 041018 (2014). 

    Download: <http://journals.aps.org/prx/pdf/10.1103/PhysRevX.4.041018>



Copyright notice
================

Copyright (c) 2014, Computational Molecular Biology Group, FU Berlin, 14195 Berlin, Germany.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

1. Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
