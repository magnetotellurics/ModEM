CSEM Version of ModEM
----------------------

.. warning::

   The CSEM version of ModEM is not fully validated and has only been used in
   some instances. We cannot guarantee correct or accurate results with it.

Please see the warning above. The CSEM version is really not yet ready for
production, but we are providing it here for future development and testing.

The CSEM version of ModEM is available on the ``CSEM`` branch of ModEM:

.. code-block:: bash

   $ git checkout CSEM

The CSEM version requires one of two submodules that are not a part of ModEM:

* Dipole1D - A part of OCCAM1DCSEM - `provided by SCRIPPS <_dipole1d>`_.
* EM1D - Provided by Rita Streich's - But not yet available for public release.

For more information on the CSEM version of ModEM see the README on the CSEM
branch: https://github.com/magnetotellurics/ModEM/blob/CSEM/README.md

.. _dipole1d: https://marineemlab.ucsd.edu/Projects/Occam/1DCSEM/


Configuring and Compiling with Dipole1D
----------------------------------------

The license associated with Dipole1D disallows us from including it in our
source code; however, the configure script (``/f90/CONFIG/configure``)
contains code to automatically retrieve it for you. The configure script has
the following options:

.. code-block:: bash

   $ cd ModEM/f90/
   $ ./CONFIG/configure 
    Usage: ./CONFIG/Configure with the following options:
    Compiler: Choose from supported compilers: [ gfortran | ifort ]
    Makefile: Provide a name for your output Makefile name.
    [Debug or Release]: Choose whether you want to compile the Debug or Release version.
    [MPI or Serial]:  Choose whether you want to compile the parallel (MPI) or serial version.
    [MF or SP or SP2]:  Choose between the Matrix Free (MF), or the Modified System of Eqs 1 (SP), or the Modified System of Eqs 2 (SP2) of the code.
    [MT or MT+CSEM]:  Compile MT or MT+CSEM. In Case of MT+CSEM, choose in the following option whether Dipole1D or EM1D or both will be used to get for the secondary field formulation.
    [Dipole1D or EM1D]:  (Optional) - Choose whether you have Dipole1D, or EM1D in the source files folder '/3D_MT/CSEM_module'

To configure with Dipole1D run:

.. code-block:: bash

   $ cd ModEM/f90
   $ ./CONFIG/configure ./CONFIG/Configure gfortran Makefile Release MPI MF MT+CSEM Dipole1D
    Dipole1D is not currently in ./3D_MT/CSEM_module/Dipole1D
    Would you like to have this script automatically download it now? [Yes/No]

You will be presented with the above prompt, type 'Yes' to automatically
download Dipole1D, then you'll be able to type ``make`` to compile:

.. code-block:: bash

   $ make

