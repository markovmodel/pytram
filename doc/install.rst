.. _ref_install:

==================
Installation guide
==================

.. toctree::
   :maxdepth: 2


There are different ways to install the pytram package.



From the repository
===================

First step: get the repository!

Go to your shell and type

.. code-block:: bash

   git clone https://github.com/markovmodel/pytram.git

Then, install the package from source.


via pip
-------

Go to the repository's root directory and type

.. code-block:: bash

   pip install .


via setup
---------

Go to the repository's root directory and type

.. code-block:: bash

   python setup.py install [--user]

To build the C-extensions in place, you can also run

.. code-block:: bash

   python setup.py build_ext --inplace



From the python package index
=============================

Go to your shell and type

.. code-block:: bash

   pip install pytram

or

.. code-block:: bash

   easy_install pytram

