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

.. note::
    
    To see the list of all CMAKE avaliable options for building ModEM, you can run:

    ``cmake . -LH``

    or if you are in ``ModEM/build``:

    ``cmake .. -LH``

Navigate into ModEM and create a new directory, ``build``. CMake is nice in that
it allows us to build out of source. Change directory into this ``build`` directory:

.. code-block:: bash

    $ cd ModEM
    $ mkdir build
    $ cd build

Now, we can run CMake:

.. code-block:: bash

    $ cmake .. -DCMAKE_Fortran_COMPILER=mpifort

The above command prints some output:

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
    -- Build files have been written to: /Users/mcurry/Projects/ModEM/build

CMake will try and locate BLAS and LAPACK automatically, if CMake has
difficulties finding them, you can specify them by using the
``LAPACK_LIBRARIES`` or the ``BLAS_LIBRARIES`` CMake command-line:

.. code-block:: bash

    $ cmake .. -DLAPACK_LIBRARIES=/path/to/lapacklib/liblapack.a -DBLAS_LIBRARIES=/path/to/lablaslib/libblas.a


After running CMake command above, and you see that the build files have been
written successfully, you'll be able to call ``make`` and build ModEM:

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


Options - Different configurations
----------------------------------

You can build different ModEM configurations by specifying different options when calling
``cmake``:

Again, once can use the following to display all the available options:

.. code-block:: bash

    $ mkdir build
    $ cd build
    $ cmake .. -LH

* 2D/3D Builds:
    * ``-DMODEM_BUILD_DIMS=<2D | 3D>``
    * Only 3D is supported.
* Forward Flavor:
    * ``-DFORWARD_FORMULATION=<MF | SP | SP2>``
* GPU Version:
    * ``-DBUILD_GPU=<CUDA | HIP>``
* MPI Flags:
    * ``-DBUILD_MPI=on/off``
* FG Flags
    * ``-DFG=on/off``
* C Timers
    * ``-DUSE_C_TIMERS=on/off``

You can specify them on the command line by doing the following during the ``cmake`` command:

.. code-block:: bash

    $ cmake .. -DFORWARD_FORMULATION=SP -DDBUILD_GPU=CUDA etc.
    $ make

.. warning::

    The ModEM CMake build can only build the 3D version of ModEM and cannot build the 2D version.

.. _findlapack_n_blas:

Can't Find LAPACK or BLAS?
---------------------------

CMake will try and locate BLAS and LAPACK automatically, if CMake has
difficulties finding them, you can specify them by using the
``LAPACK_LIBRARIES`` or the ``BLAS_LIBRARIES`` CMake command-line:

.. code-block:: bash

    $ cmake .. -DLAPACK_LIBRARIES=/path/to/lapacklib/liblapack.a -DBLAS_LIBRARIES=/path/to/lablaslib/libblas.a


Specifying Intel Instruction Set
---------------------------------

If you are building with the Intel compiler, you may want to specify the feature
set, to do this, you can specify additional compler flags using ``-DCMAKE_Fortran_FLAGS```:

.. code-block:: bash

    $ cmake .. -DCMAKE_Fortran_FLAGS="-xHost"

For the Ifort compiler, ``-xHost`` will tell the compiler to generate the
highest instruction set available for your host process. As well, you can
also explicitly specify your feature set:

.. code-block:: bash

    $ cmake .. -DCMAKE_Fortran_FLAGS="-msse4.2"



Building ModEM Spherical Version
==================================

.. warning::

    The Spherical version of ModEM is experimental and will not work with
    traditional MT source setup. MT forward modeling in spherical coordinates
    relies on an external 1D forward model and the secondary field formulation
    option in ModEM to be physically correct.

    Using spherical coordinates with the tradtional MT source setup has not
    been shown to be physically accurate.


To build the experimental spherical version of ModEM you can specify the
``BUILD_SPHERICAL`` argument in CMake:

.. code-block:: bash

   $ cmake .. -DBUILD_SPHERICAL=on


After running the above, you will presented with the same warning printed
above, in order to fully build the experimental version, you must also pass the
`FORCE_SPHERICAL` argument as well:

.. code-block:: bash

   $ cmake .. -DBUILD_SPHERICAL=on -DFORCE_SPHERICAL=on
   $ make


Building GPU capable ModEM with CMake
======================================

The nice thing about CMake over the old configure script is that CMake can
compile languages other than just Fortran. This means that it can help us
compile the GPU code quite easily. Here are the instructions on how to compile 
the CUDA code.


Compiling ModEM GPU Code with CUDA
-----------------------------------

You will need to have the `CUDA toolkit installed <cuda-toolkit_>`_, so far,
ModEM has been confirmed to build and work with cuda11.3/toolkit. You may try
other version, but your results may vary. Of course, you will also need to have
a Nvidia capable GPU installed as well.

The toolkit will provide a number of CUDA libraries that ModEM requires: cublas,
cudart, cusparse and NVML (Nvidida Management Library). Thankfully, CMake will
help us with linking this, you just need to make sure that ``CUDA_PATH``
environment is set to the toolkit installation location and/or ``nvcc`` is
available on the command line.

To build ModEM with GPU CUDA we can run the following command:

.. code-block:: bash

    $ mkdir build; cd build
    $ cmake .. -DCMAKE_Fortran_COMPILER=mpifort -DBUILD_GPU=CUDA

Again, you may need to specify LAPACK or BLAS libraries, see
:ref:`findlapack_n_blas`. After building, you should see confirmation that the
CUDA compiler (``nvcc``) has been found and that the CUDAToolkit was found:

.. code-block:: bash

    -- The CUDA compiler identification is NVIDIA 11.3.58
    -- Detecting CUDA compiler ABI info
    -- Detecting CUDA compiler ABI info - done
    -- Check for working CUDA compiler: /cm/shared/apps/cuda11.3/toolkit/11.3.0/bin/nvcc - skipped
    -- Detecting CUDA compile features
    -- Detecting CUDA compile features - done
    ...
    -- Found CUDAToolkit: /cm/shared/apps/cuda11.3/toolkit/11.3.0/include (found version "11.3.58")
    ...
    -- Configuring done
    -- Generating done
    -- Build files have been written to: /home/mcurry/Projects/ModEM/build

If CMake returned success, you will be able to call ``make`` and build the CUDA compatible ModEM build:

.. code-block:: bash

    $ make
    [  1%] Building Fortran object f90/CMakeFiles/Mod3DMT_SP2.dir/UTILS/math_constants.f90.o
    [  3%] Building Fortran object f90/CMakeFiles/Mod3DMT_SP2.dir/UTILS/utilities.f90.o
    [  5%] Building Fortran object f90/CMakeFiles/Mod3DMT_SP2.dir/3D_MT/DICT/dataTypes.f90.o
    ... # Ommitting build output...
    [ 96%] Building Fortran object f90/CMakeFiles/Mod3DMT_SP2.dir/Mod3DMT.f90.o
    [ 98%] Building CUDA object f90/CMakeFiles/Mod3DMT_SP2.dir/3D_MT/FWD_SP2/kernel_c.cu.o
    [100%] Linking Fortran executable Mod3DMT_SP2
    [100%] Built target Mod3DMT_SP2

Specifying CUDA GPU Compute Capability
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When building with CUDA GPU code you may want to specify what `CUDA GPU Compute Capability <cuda-compute-capability_>`_ 
you build with. It is probably best to target the highest compute capability
that your GPU supports.

By default, NVCC will target the lowest compute capability possible. To find
your supported compute capability version, you can use the ``nvidida-smi``
command to query your installed GPU:

.. code-block::

    $ nvidia-smi --query-gpu=name,compute_cap --format=csv
    name, compute_cap
    Tesla V100-SXM2-16GB, 7.0

You can specify your desired compute capability by specifying the
``CMAKE_CUDA_ARCHITECTURES`` CMake command line option:

.. code-block:: bash

    $ cmake .. -DCMAKE_Fortran_COMPILER=mpifort -DBUILD_GPU=CUDA -DCMAKE_CUDA_ARCHITECTURES="35;50;72"

This will generate code for Cuda architectures 35, 50 and 72. You can also tell
CMake to use whatever capability your GPU can handle with the ``native`` keyword:

.. code-block:: bash

    $ cmake .. -DCMAKE_Fortran_COMPILER=mpifort -DBUILD_GPU=CUDA -DCMAKE_CUDA_ARCHITECTURES="native"

You can also specify ``all`` or ``all-major``, which will generate code for all
or all major version. For more information on available options see:
https://cmake.org/cmake/help/latest/prop_tgt/CUDA_ARCHITECTURES.html#prop_tgt:CUDA_ARCHITECTURES.


.. _cuda-compute-capability: https://developer.nvidia.com/cuda/gpus

Running CUDA ModEM 
~~~~~~~~~~~~~~~~~~~

When running with the GPU code, you should see output that confirms your GPU has been found:

.. code-block::

    node[001]:    Waiting for a message from Master
    node[001]:    total time elapsed (sec):     0.000000
    node[001]:    MPI TASK [                       FORWARD] received from     0
    node[001]:  Computing the BC internally  for period    1      & mode #  1
    node[001]:  Solving the FWD problem for period    1:  1.163636E+01 secs & mode #  1
    # Dev Status  : GPU-mem = 1.916719 %, PCI 98: 0
    # Dev Selected: 0. Tesla V100-SXM2-16GB

However, it should be noted that ModEM does not fail if it can't detect your
GPU, it just falls back to the CPU version of the code:

.. code-block::

    node[001]:    Waiting for a message from Master
    node[001]:    total time elapsed (sec):     0.000000
    node[001]:    MPI TASK [                       FORWARD] received from     0
    node[001]:  Computing the BC internally  for period    1      & mode #  1
    node[001]:  Solving the FWD problem for period    1:  1.163636E+01 secs & mode #  1
    [WARNING] could not find a valid GPU...
    [WARNING] Fall back to CPU version of BICG

If this is the case, you should ensure that the ``nvidida-smi`` command line tool can find your GPU:

.. code-block:: bash

    $ nvidida-smi
    +-----------------------------------------------------------------------------------------+
    | NVIDIA-SMI 550.54.15              Driver Version: 550.54.15      CUDA Version: 12.4     |
    |-----------------------------------------+------------------------+----------------------+
    | GPU  Name                 Persistence-M | Bus-Id          Disp.A | Volatile Uncorr. ECC |
    | Fan  Temp   Perf          Pwr:Usage/Cap |           Memory-Usage | GPU-Util  Compute M. |
    |                                         |                        |               MIG M. |
    |=========================================+========================+======================|
    |   0  Tesla V100-SXM2-16GB           On  |   00000000:89:00.0 Off |                    0 |
    | N/A   28C    P0             41W /  300W |       0MiB /  16384MiB |      0%      Default |
    |                                         |                        |                  N/A |
    +-----------------------------------------+------------------------+----------------------+

    +-----------------------------------------------------------------------------------------+
    | Processes:                                                                              |
    |  GPU   GI   CI        PID   Type   Process name                              GPU Memory |
    |        ID   ID                                                               Usage      |
    |=========================================================================================|
    |  No running processes found                                                             |
    +-----------------------------------------------------------------------------------------+


.. _cuda-toolkit: https://developer.nvidia.com/cuda/toolkit

Compiling ModEM GPU Code with HIP
----------------------------------

Compiling the HIP code with CMAKE is currently a WIP.


