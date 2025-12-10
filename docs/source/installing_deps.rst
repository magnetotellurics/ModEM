
.. _installing-deps:

Installing ModEM Dependencies
=============================

Installing MPI
--------------

There are two widely used versions of MPI which we recommend using with ModEM:

* MPICH - https://www.mpich.org/
* OpenMPI - https://www.open-mpi.org/

Both have been tested throughly with ModEM. Both offer benifits over ther 
another, but such disccusion is beyond the scope of this document. Installing
both is reletivley similar, but for brevity we will only give instructions on
installing MPICH. 

Linux
^^^^^^

For Linux distributions you can install it most often via your package panager
(if you have sudo powers). i.e.:

.. code-block:: bash

    # Debian/Ubuntu
    $ apt install openmpi-bin libopenmpi-dev 

    # Fedora/RHEL/Centos
    $ yum install mpich

.. note::

    Installing with a Linux package manager requires root (``sudo``) powers,
    thus the above commands will not work on a shared system where you do not
    have powers to use ``sudo``.

If you do not have sudo/root powers, then you can follow the instructions to
compile MPICH by hand in :ref:`mpich-by-hand`.

MacOS
^^^^^^


For MacOSX, you can also install using Homebrew:

.. code-block:: bash

    $ brew install mpich

See: https://formulae.brew.sh/formula/mpich

.. _mpich-by-hand: 

Compiling mpich from source 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. tip::

    More information on install MPICH (or OpenMPI) can be found in 
    thier README or INSTALL files of their souce code.

It is relativly easy to install MPICH (and OpenMPI) by hand if you
need to do so. To do so, you will only need to a C, C++, and Fortran 
compiler avaiable on your system.

First, download the zip/tar file from the MPICH website and then unzip it. 

.. code-block:: bash

    $ wget https://www.mpich.org/static/downloads/5.0.0b1/mpich-5.0.0b1.tar.gz
    $ tar -xzvf mpich-5.0.0b1.tar.gz

Then, create a new, empty, directory inside your home directory, this is where
we will install our files too.

.. code-block:: bash

    $ cd ~/
    $ mkdir installs

We will use the ``installs`` directory to install MPI, and the LAPACK and
BLAS (as described below). 

Change directory into the MPICH direcotry:

.. code-block:: bash

    $ cd mpich-5.0.0.b1

Then, we can run the ``configure`` inside the mpich directory. Here, we will specify
our C compiler and Fortran compiler (although they are often detected by default):

.. code-block:: bash

    $ CC=gcc FC=gfortran ./configure --prefix=/home/<USERNAME>/installs

In the line above, CC sets the C compile we want to use (to gcc) and FC sets the
Fortran compiler we want to use (gfortran). You can adjust these as you see fit. 

Configure will take sometime to determine how to build your system, and at the
end, you'll get a message saying that 'MPICH is configured..'. Once it's configured,
we can call ``make`` to compiler it:

.. code-block:: bash

    $ # Make with 4 threads
    $ make -j4

Here we've used ``-j`` to specify the number of threads to use when compiling,
you can use more, but if you are on a shared system you should respect other
users and use an appropriate amount of threads (4 is generally good).

Once it's compiled, you can run ``make install``, which will install MPICH into
the the ``/home/<USERNAME>/installs``

.. code-block:: bash

    $ make -j4 install


Inside ``/home/<USERNAME>/installs`` there will be four new directories:

* ``bin`` - Contains executables such as ``mpifort``, ``mpicc``, ``mpiexec``, ``mpirun`` etc.
* ``etc`` - Contains MPI configuration files
* ``include`` - Contains MPI .h and .mod include files
* ``lib`` - MPI library files
* ``share`` - MPI Documents and manuals ``man mpifort``.

To use our install, we will need to run the executables that are in the
``/home/<USERNAME>/installs/bin`` directory. To do so we can add them to our PATH
variable:

.. code-block:: bash

    $ export PATH=/home/<USERNAME>/install/bin:$PATH

Now we can run ``mpifort`` to compile and ``mpiexec`` to launch MPI applications.

Installing LAPACK and LBLAS
----------------------------

Installing BLAS
^^^^^^^^^^^^^^^^

While you can install BLAS manually, LAPACK *is* packaged with a copy of BLAS, so
you can skip straight to installing LAPACK: See :ref:`installing-lapack`.


.. _installing-lapack:

Installing LAPACK
^^^^^^^^^^^^^^^^^^

.. tip::

    More information for installing LAPACK can be found in its README.md file.

Installing LAPACK is very easy. First, download a copy of LAPACK from:
https://www.netlib.org/lapack/ and extract it:

.. code-block:: bash

    $ wget https://github.com/Reference-LAPACK/lapack/archive/refs/tags/v3.12.1.tar.gz
    $ tar -xzvf v3.12.1.tar.gz
    $ cd lapack-3.12.1

For LAPACK, we just need to create a ``make.inc`` file. We simply can use the
provided example ``make.inc.example``:

.. code-block:: bash

    $ cp make.inc.example make.inc

Open ``make.inc`` and insure it has correct settings (it should be set for gcc/gfortran),
then call make:

.. code-block:: bash

    $ make -j4

This will produce three ``*.a`` files:

* liblapack.a
* librefblas.a
* libtmglib.a

LAPACK does not have an install like MPI, but we can just copy, link or move
them into the installation we made in :ref:`mpich-by-hand`
(``/home/<USERNAME>/installs/lib``). Specificlaly, we will need to move/link them
into the ``lib`` folder of that directory.

We will also need to change the name of ``librefblas.a`` to ``libblas.a``, as ModEM
tries to link to ``libblas.a``.

.. code-block:: bash

    $ ln -s liblapack.a /home/<USERNAME>/install/lib/liblapack.a
    $ ln -s librefblas.a /home/<USERNAME>/install/lib/libblas.a


Compiling ModEM with LAPACK/LABLAS
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Now, when we compile ModEM, we can link to these libraries. To do so, we will
need to specify the directory ``lib`` directory above, we can do this in the
``LDFLAGS`` when we run configure:

.. code-block:: bash

    $ cd ModEM/f90
    $ LDFLAGS="-L/home/<USERNAME>/installs/lib" ./CONFIG/configure Makefile gfortran

This will cause your installed LAPACK and BLAS to be linked to ModEM!
