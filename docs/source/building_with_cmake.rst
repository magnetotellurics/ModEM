Building ModEM With CMake
=========================

ModEM now has the ability to be compiled with CMake.


Install CMake
--------------

You can download CMake from their website: https://cmake.org/download/. 

For MacOSX you can install CMake via homebrew

.. code-block:: bash

    $ brew install cmake

On Linux, you should be able to install CMake in the same way using the package
manager. I'm sure google can help with this...

Note About CMake Caching
------------------------

Confusing to many new users, CMake is very eager to cache the arguments that you
send it, thus, when you want to use different options, you will need to remove
these cached variables, which is in ``CMakeCache.txt``.

If you are having strange troubles switching between try removing
``CMakeCache.txt`` and running ``cmake`` again.


Building ModEM with CMake
--------------------------

Navigate into ModEM and create a new directory, ``build``. CMake is nice in that
it allows us to build out of source. Change directory into this ``build`` directory:


.. code-block:: bash

    $ cd ModEM
    $ mkdir build
    $ cd build

Now, we can run CMake:

.. code-block:: bash

    $ cmake .. -DCMAKE_Fortran_COMPILER=mpifort
    -- The Fortran compiler identification is LLVMFlang 20.1.8
    -- Detecting Fortran compiler ABI info
    -- Detecting Fortran compiler ABI info - done
    -- Check for working Fortran compiler: /opt/homebrew/bin/flang - skipped
    -- Building BUILD_MF
    -- Adding sources too Mod3DMT
    -- Building ModEM ExE: Mod3DMT
    -- Looking for Fortran sgemm
    -- Looking for Fortran sgemm - not found
    -- Looking for Fortran dgemm
    -- Looking for Fortran dgemm - not found
    -- Looking for Fortran dgemm
    -- Looking for Fortran dgemm - not found
    -- Looking for Fortran sgemm
    -- Looking for Fortran sgemm - found
    -- Found BLAS: /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib/libblas.tbd
    --  WE FOUND BLAS
    -- Looking for Fortran cheev
    -- Looking for Fortran cheev - not found
    -- Looking for Fortran cheev
    -- Looking for Fortran cheev - not found
    -- Looking for Fortran cheev
    -- Looking for Fortran cheev - found
    -- Found LAPACK: /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib/liblapack.tbd;/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib/libblas.tbd
    --  WE FOUND LAPACK
    -- Configuring done (2.1s)
    -- Generating done (0.0s)
    -- Build files have been written to: /Users/mcurry/Projects/ModEM-Model/build

There may be some additonal warnings that I have removed for breviity. 

CMake will try and locate BLAS and LAPACK automatically, if you have installed
them manually, I belive you can specify their location by doing:

.. code-block:: bash

    $ # But I still need to test this...
    $ cmake .. -DCMAKE_EXE_LINKER_FLAGS="-L /path/to/lblas"


After running CMake above, and you see the last line: ``-- Build files have been
written to: /Users/mcurry/Projects/ModEM-Model/build`` you will be ready to call
``make`` and build ModEM:

.. code-block:: bash

    $ make
    [  1%] Building Fortran object f90/CMakeFiles/Mod3DMT.dir/UTILS/math_constants.f90.o
    [  3%] Building Fortran object f90/CMakeFiles/Mod3DMT.dir/UTILS/utilities.f90.o
    [  5%] Building Fortran object f90/CMakeFiles/Mod3DMT.dir/3D_MT/DICT/dataTypes.f90.o
    [  7%] Building Fortran object f90/CMakeFiles/Mod3DMT.dir/3D_MT/DICT/receivers.f90.o
    [  9%] Building Fortran object f90/CMakeFiles/Mod3DMT.dir/3D_MT/DICT/transmitters.f90.o
    [ 11%] Building Fortran object f90/CMakeFiles/Mod3DMT.dir/UTILS/elements.f90.o
    [ 13%] Building Fortran object f90/CMakeFiles/Mod3DMT.dir/3D_MT/GridDef.f90.o
    [ 15%] Building Fortran object f90/CMakeFiles/Mod3DMT.dir/FIELDS/FiniteDiff3D/sg_scalar.f90.o
    [ 16%] Building Fortran object f90/CMakeFiles/Mod3DMT.dir/FIELDS/FiniteDiff3D/sg_vector.f90.o
    [ 18%] Building Fortran object f90/CMakeFiles/Mod3DMT.dir/FIELDS/FiniteDiff3D/sg_spherical.f90.o
    [ 20%] Building Fortran object f90/CMakeFiles/Mod3DMT.dir/3D_MT/GridCalc.f90.o
    [ 22%] Building Fortran object f90/CMakeFiles/Mod3DMT.dir/UTILS/file_units.f90.o
    [ 24%] Building Fortran object f90/CMakeFiles/Mod3DMT.dir/FIELDS/FiniteDiff3D/sg_sparse_vector.f90.o
    [ 26%] Building Fortran object f90/CMakeFiles/Mod3DMT.dir/3D_MT/modelParam/ModelSpace.f90.o
    [ 28%] Building Fortran object f90/CMakeFiles/Mod3DMT.dir/3D_MT/EMfieldInterp.f90.o
    [ 30%] Building Fortran object f90/CMakeFiles/Mod3DMT.dir/UTILS/fields_orientation.f90.o
    [ 32%] Building Fortran object f90/CMakeFiles/Mod3DMT.dir/FIELDS/FiniteDiff3D/sg_boundary.f90.o
    [ 33%] Building Fortran object f90/CMakeFiles/Mod3DMT.dir/3D_MT/SolnSpace.f90.o
    [ 35%] Building Fortran object f90/CMakeFiles/Mod3DMT.dir/3D_MT/DataFunc.f90.o
    [ 37%] Building Fortran object f90/CMakeFiles/Mod3DMT.dir/SENS/DataSpace.f90.o
    [ 39%] Building Fortran object f90/CMakeFiles/Mod3DMT.dir/3D_MT/DataIO.f90.o
    [ 41%] Building Fortran object f90/CMakeFiles/Mod3DMT.dir/3D_MT/FWD/Mod2d/WSfwd2Dpar.f90.o
    [ 43%] Building Fortran object f90/CMakeFiles/Mod3DMT.dir/3D_MT/FWD/Mod2d/WSutils.f90.o
    [ 45%] Building Fortran object f90/CMakeFiles/Mod3DMT.dir/3D_MT/FWD/Mod2d/WSfwd1Dmod.f90.o
    [ 47%] Building Fortran object f90/CMakeFiles/Mod3DMT.dir/3D_MT/FWD/Mod2d/WSfwd2Dmod.f90.o
    [ 49%] Building Fortran object f90/CMakeFiles/Mod3DMT.dir/3D_MT/FWD/Mod2d/FwdTEmod.f90.o
    [ 50%] Building Fortran object f90/CMakeFiles/Mod3DMT.dir/3D_MT/FWD/boundary_ws.f90.o
    [ 52%] Building Fortran object f90/CMakeFiles/Mod3DMT.dir/3D_MT/FWD/nestedEM.f90.o
    [ 54%] Building Fortran object f90/CMakeFiles/Mod3DMT.dir/FIELDS/FiniteDiff3D/sg_diff_oper.f90.o
    [ 56%] Building Fortran object f90/CMakeFiles/Mod3DMT.dir/3D_MT/FWD/modelOperator3D.f90.o
    [ 58%] Building Fortran object f90/CMakeFiles/Mod3DMT.dir/3D_MT/FWD/solver.f90.o
    [ 60%] Building Fortran object f90/CMakeFiles/Mod3DMT.dir/3D_MT/FWD/EMsolve3d.f90.o
    [ 62%] Building Fortran object f90/CMakeFiles/Mod3DMT.dir/3D_MT/FWD/Mod2d/FwdTMmod.f90.o
    [ 64%] Building Fortran object f90/CMakeFiles/Mod3DMT.dir/3D_MT/ForwardSolver.f90.o
    [ 66%] Building Fortran object f90/CMakeFiles/Mod3DMT.dir/3D_MT/ioMod/ioAscii.f90.o
    [ 67%] Building Fortran object f90/CMakeFiles/Mod3DMT.dir/SENS/SensMatrix.f90.o
    [ 69%] Building Fortran object f90/CMakeFiles/Mod3DMT.dir/UserCtrl.f90.o
    [ 71%] Building Fortran object f90/CMakeFiles/Mod3DMT.dir/3D_MT/Main.f90.o
    [ 73%] Building Fortran object f90/CMakeFiles/Mod3DMT.dir/3D_MT/SolverSens.f90.o
    [ 75%] Building Fortran object f90/CMakeFiles/Mod3DMT.dir/3D_MT/Sub_MPI.f90.o
    [ 77%] Building Fortran object f90/CMakeFiles/Mod3DMT.dir/SENS/DataSens.f90.o
    [ 79%] Building Fortran object f90/CMakeFiles/Mod3DMT.dir/SENS/SensComp.f90.o
    [ 81%] Building Fortran object f90/CMakeFiles/Mod3DMT.dir/INV/DCG.f90.o
    [ 83%] Building Fortran object f90/CMakeFiles/Mod3DMT.dir/INV/INVcore.f90.o
    [ 84%] Building Fortran object f90/CMakeFiles/Mod3DMT.dir/INV/LBFGS.f90.o
    [ 86%] Building Fortran object f90/CMakeFiles/Mod3DMT.dir/UTILS/ModEM_timers.f90.o
    [ 88%] Building Fortran object f90/CMakeFiles/Mod3DMT.dir/INV/NLCG.f90.o
    [ 90%] Building Fortran object f90/CMakeFiles/Mod3DMT.dir/MPI/Declaration_MPI.f90.o
    [ 92%] Building Fortran object f90/CMakeFiles/Mod3DMT.dir/MPI/Main_MPI.f90.o
    [ 94%] Building Fortran object f90/CMakeFiles/Mod3DMT.dir/SENS/SymmetryTest.f90.o
    [ 96%] Building Fortran object f90/CMakeFiles/Mod3DMT.dir/UTILS/polpak.f90.o
    [ 98%] Building Fortran object f90/CMakeFiles/Mod3DMT.dir/Mod3DMT.f90.o
    [100%] Linking Fortran executable Mod3DMT
    [100%] Built target Mod3DMT

The ``Mod3DMT`` will then be placed in the f90 file inside your build directory,
which you will be able to run as normal.

Showing CMake Options
----------------------

You can show the options aviable in the CMake build by running:

.. code-block:: bash

    $ cmake -LH

That will show you the list of options that you can build with ModEM, as well as
their help message. It will also show you some additional help options as well.

Options - Different configurations
----------------------------------

You can build different ModEM configurations by specifying different options when calling
``cmake``:

.. warning::

    These flags currently do not do anything!

* 2D/3D Builds:
    * ``-DMODEM_BUILD_DIMS=<2D | 3D>``
* Forward Flavor:
    * ``-DFORWARD_FLAVOR=<MF | SP | SP2>``
* MPI Flags:
    * ``-DBUILD_MPI=on/off``
* FG Flags
    * ``-DFG=on/off``
* Flags Currently not added: 
    * ``-DFG=on/off``
    * ``-DHIP=on/off``
    * ``-DCUDA=on/off``

You can specify them on the command line by doing the following during the ``cmake`` command:

.. code-block:: bash

    $ cmake .. -DFORWARD_FLAVOR=SP
    $ make
