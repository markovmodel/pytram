.. _ref_install:

==================
Installation guide
==================

.. toctree::
   :maxdepth: 2


There are different ways to install the pytram package.

From the python package index
=============================

Go to your shell and type

.. code-block:: bash

   pip install pytram

or

.. code-block:: bash

   easy_install pytram

Possible pitfalls can be the 'Cython' dependency. In case the installation fails,
try to install Cython via the package index as well, by typing into your shell:

.. code-block:: bash

   pip install cython


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

From the repository
===================

First step: get the repository!

Go to your shell and type

.. code-block:: bash

   git clone https://github.com/markovmodel/pytram.git

Then, install the package from source.

If you expereince any other difficulties with the installation, please email us using the mailing list pytram@lists.fu-berlin.de.
Ideally provide us with some basic information about your operating system/python distribution, so that we can try and reproduce your issues. 



