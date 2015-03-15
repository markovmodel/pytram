.. _ref_user:

==========
User guide
==========

.. toctree::
   :maxdepth: 2


Getting started
===============

The pytram package comes with two examples for dTRAM and xTRAM in the form of ipython notebooks. In order to run run them, go the root directory of the pytram repositroy and type

.. code-block:: bash

   ipython notebook dTRAM_example.ipynb

in your shell for dTRAM or 

.. code-block:: bash

   ipython notebook xTRAM_example.ipynb

for xTRAM, respectively. These examples will illustrate the basic usage of the API functions and the necessary data preparation.

Another way is a file-based approach which allows you to store simulation results for later analysis. To this aim, we have specified the required file format in the next section and implemented a reader function to import such data into pytram. Using file-based input data allows you to run pytram

  * from the API using a ``Reader`` and the ``TRAMData`` converter or
  * using run scripts directly from the shell.

The latter option requires the least effort, however, also offers only limited flexibility.


File format
===========

The standard file format assumes text files with the following layout. ::

   # This is a comment line, you can several of those.
   # The next lines indicates the meaning of the columns,
   # [M] denotes Markov state indices (starting from zero),
   # [T] denotes thermodynamic state indices (starting from zero),
   # and [b_K] denotes the reduced bias energies b_K/kT
   # [M]  [T]  [b_0]  [b_1]  ... 
      0    0    3.1   18.2
      1    0    3.2   18.3
      2    0    4.8   19.9
      3    0    7.4   22.5
      .    .      .      .
      .    .      .      .
      .    .      .      .

The minimal layout only requires the ``[M]`` and ``[T]`` columns and can only be used for dTRAM. These two columns contain the sequences of the Markov and generating thermodynamic states. For example, the entry ``3  0`` denotes that the actual sample corresponds to the Markov state ``3`` and was generated at thermodynamic state ``0``.

**Important note**: in order to run dTRAM successfully, you need an additional ``b_K_i`` file as explained in the dTRAM section.

The other TRAM estimators require at least one energy column. For this, we distinguish two different cases:
temperature as the only thermodynamic variable, and all other thermodynamic conditions, i.e., different Hamiltonians, umbrella potentials, ...


Temperature as only thermodynamic variable
------------------------------------------

In this case, you need the ``[M]`` and ``[T]`` columns, and one energy column ``[b_K]``; this column contains the reduced energy sequence. The energy is reduced according to the generating thermodynamic state. For example, the entry ``2  5  20.5`` denotes that the actual sample corresponds to the Markov state ``2``, was generated at temperature ``kT_5``, and the corresponding energy was reduced with ``kT_5``.

**Important note**: for temperature-dependent simulations, you need an additional single column ``kT`` file wich indicates all generating temperatures times the Boltzmann constant (consistent with your energy units). Note that the order of ``kT`` values must be constistent with the numbering of the thermodynamic states.


Hamiltonian replica exchange, umbrella sampling, etc
----------------------------------------------------

This is the most general application. Here, each sample must be evaluated at all thermodynamic states which means that you need as many energy columns as you have thermodynamic states. For example, the line ``2  1  3.0  2.9  1.0  0.3`` indicates that the actual sample corresponds to the Markov state ``2``, has been generated at thermodynamic state ``1``, the reduced energy is

  * ``3.0 kT`` at thermodynamic state ``0``,
  * ``2.9 kT`` at thermodynamic state ``1``,
  * ``1.0 kT`` at thermodynamic state ``2``, and
  * ``0.3 kT`` at thermodynamic state ``3``.

This example also requires you to have exactly four thermodynamic states.


Running dTRAM
=============

from files
----------

Assume that we have two files ``file_1.dat`` and ``file_2.dat`` with simulation data. In addition to that, the dTRAM method requires the user to specify the reduced bias energies of all Markov states in each of the thermodynamic states. The corresponding file format is given by ::

   # we store the reduced bias energies b_K(x)/kT
   # at the discrete states x_i
   # [b_0]  [b_1]  ... 
      0.0    4.0
      0.0    0.0
      0.0    8.0

In this example, we have three Markov states which are evaluated at two different thermodynamic states.

Using the API, we can run dTRAM via the following code:

.. code-block:: python

   # import the Reader, TRAMData and the dtram API function
   from pytram import Reader, TRAMData, dtram

   # specify your input data files
   files = [ 
         'path/to/file_1.dat',
         'path/to/file_2.dat'
      ]
   b_K_i_file = 'path/to/b_K_i_file.dat'

   # import the files using the Reader
   reader = Reader( files, b_K_i_file=b_K_i_file, verbose=True )

   # convert the input data using TRAMData
   data = TRAMData( reader.trajs, b_K_i=reader.b_K_i )

   # run dTRAM using the API function
   dtram_obj = dtram( data, maxiter=1000, ftol=1.0E-10, verbose=True )

   # show unbiased stationary distribution
   print dtram_obj.pi_i

   # get transition matrix for thermodynamic state 0
   T_0 = dtram_obj.estimate_transtition_matrix( 0 )
   print T_0

   # show thermodynamic free energies
   print dtram_obj.f_K

Optionally, we can use the dTRAM run script from shell. Simply type

.. code-block:: bash

   dtram.py --b_K_i_file=path/to/b_K_i_file.dat --verbose path/to/file_1.dat path/to/file_2.dat

and enjoy the show. Note that the run script will not work if the pytram package has not been properly installed!

You can run

.. code-block:: bash

   dtram.py --help

for additional information on available parameters.


from seqential data
-------------------

The data preparation and the API usage is shown in the ipython example.





