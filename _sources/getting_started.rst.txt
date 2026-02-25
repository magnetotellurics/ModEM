Getting Started
================

Directory Structure
--------------------

The following files and sub-directories will be found:

* **COPYRIGHT** - Please get familiar with the Copyright before using this code!
* **README.md** - Explains how to obtain, install and update the code.
* **doc/** - Provides additional documentation, including this user guide.
* **f90/** - The code base, makefiles and configuration scripts.

Compiling ModEM
---------------

The build system for ModEM uses a ModEM configure script and the ``fmkmf.pl``
Perl script:

* configure script (``./f90/config/configure``) - Calls ``fmkmf.pl`` with
  various arguments
* ``fmkmf.pl`` - Called by Configure scripts and automatically generates a
  ModEM makefile.

Once you have run configure, you'll be able to run ``make`` on the generated
Makefile.

In the ``f90`` directory there are some prototype Makefiles which can be used
to compile ModEM without having to run the configure script; however, these
Makefiles only create the Matrix-Free version of ModEM and you may want to
compile the other two version SP/SP2. To do this you will need to run the
configure script, which is explained below.

ModEM Dependencies
------------------

ModEM requires two libraries, the BLAS Library and the LAPACK library and
optionally the MPI library.

Often BLAS and LAPACK are already available on your system and are linked
automatically, but occasionally you might need to install or link them
yourself. For instructions on building ModEM dependencies see
:ref:`installing-deps`.

Configure Scripts
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

There are a number of configure scripts inside ``./f90/CONFIG`` in the form of
``Configure.***``, ``Configure.SP.***``, ``Configure.3D_MT.**``, etc. These
configure files are considered deprecated, and the ``configure`` script inside
``f90/CONFIG`` should be used instead; however, the deprecated configure
scripts are left for backwards compatibility.

Creating Makefiles from configure 
----------------------------------

.. note::

    By default, ``configure``, with no optional arguments will generate an MPI
    Makefile with the SP2 forward solver.

To create a Makefile, run the ``configure`` script found in the ``f90`` with
arguments that match your system and the version of ModEM you want to use.

Currently, ModEM has three different forwards 'flavors':

* MF - Matrix Free - Original/classic ModEM
* SP - Hao's Sparse Matrix
* SP2 -  (DEFAULT) - Hao's Second Sparse Matrix - Fastest and has Hao's
    fine-grained capability

.. code-block:: bash

    $ cd f90/
    $ # By default ./CONFIG/configure creates a SP2 Makefile with MPI
    $ ./CONFIG/configure Makefile gfortran

The Configure scripts has the following usage::

    $ Usage: ./CONFIG/Configure [-h] [-u] [-g release|debug] [-m mpi|serial] [-l solver] [-f] makefile-name compiler
        -h Display full help message
        -u Display this usage message
        -g <release | debug> Specify whether to build with debug symbols or not
            release - Default - Create Makefile with full optimizations (i.e. -O3)
            debug   - Create Makefile with debug symbols (i.e. -g) and no optimizations
        -m <mpi* | serial> Specify MPI or serial
            mpi - Default - Builds with MPI
            serial - Builds Serially
        -l <MF | SP | SP2> Specify solver
            MF - Builds with MF
            SP - Builds with SP
            SP2 - Default - Builds with SP2
        -f Build with Fined-Grained -DFG on SP2
            Only able to build with SP2
        -s Build spherical version of ModEM

Running configure will create a Makefile that you specify in the first
argument. In the above case the Makefile will be named ``Makefile``. We can
then call ``make`` on that ``Makefile``

.. code-block:: bash

    $ make

Examples Using Configure
^^^^^^^^^^^^^^^^^^^^^^^^^

A few examples using Configure

.. code-block:: bash

    $ ./CONFIG/configure Makefile gfortran
    $ # Generates a Makefile with GFortran options using MPI

    $ ./CONFIG/configure -d debug Makefile ifort
    $ # Generates an MPI Makefile with ifort debugging options

    $ ./CONFIG/configure -m serial -l MF Makefile gfortran
    $ # Generates a serial Makefile with GFortran with the MF solver

    $ FC=ftn LDFLAGS=-L/home/username/install/lib ./CONFIG/configure Makefile gfortran
    $ # Generates an SP2, MPI capabile makefile, adding LDFLAGS to the link step and
    $ # specifies the compiler to be `ftn`.

    $ FFLAGS='-mSSE2' ./CONFIG/configure Makefile ifort
    $ # Create a SP2, MPI capabile makefile, with ifort/intel compiler options and add
    $ # the '-mSSE2' to all compiling steps.

    $ LDFLAGS="-L/path/to/install/lib/" ./CONFIG/configure Makefile gfortran
    $ # Generates a Makefile with GFortran options using MPI and passes
    $ # the location of /path/to/install/lib to the linker step

Configuring on Cray Systems
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you are running on a Cray system, you will need to use Cray's Fortran
compiler wrapper: ``ftn``. You can easily specify the name of the compiler
by setting the ``FC`` environment variable:

.. code-block:: bash

    $ FC=ftn ./CONFIG/configure Makefile gfortran
    $ # Or for intel:
    $ FC=ftn ./CONFIG/configure Makefile ifort


Ensure BLAS and LAPACK are linked
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you get errors that you are not able to find ``-lblas`` or ``-llapack``, but
have them installed, you can use the ``LDFLAGS`` environment variable to tell
the linker where to find them. You can specify this when you run configure:

.. code-block:: bash

   $ LDFLAGS="-L/path/to/install/lib/" ./CONFIG/configure Makefile gfortran

If you have already created your makefile, you can also add the ``-L``
specification to the ``LIBS_PATH`` of the Macro-Def section of your Makefile:

.. code-block:: Make

    # ------------------Macro-Defs---------------------
    include Makefile.local
    OBJDIR = ./objs/3D_MT/csemBuild
    F90 = mpifort
    FFLAGS = -cpp -DMPI -DCSEM_Dipole1D -O3 -ffree-line-length-none -dI -fallow-argument-mismatch
    MPIFLAGS = -cpp -DMPI -DCSEM_Dipole1D
    MODULE = -J $(OBJDIR)
    LIBS_PATH = -L/path/to/install/lib
    LIBS = -lblas -llapack


Compiling
----------

Once you have generated your own makefile, or decided to use one of the
defaults, you can compile ModEM by renaming that makefile to Makefile and then
run make

.. code-block:: bash

    $ cp Makefile.3D.MF.gnu Makefile
    $ make clean # Not necessary from a fresh clone
    $ make

Compiling with MPI
^^^^^^^^^^^^^^^^^^

If you have MPI installed, and already have a serial makefile generated you can
easily compile the MPI version of ModEM by altering the makefile. This is
sometimes easier then creating a new makefile from the configuration scripts.

To compile an MPI version, change the ``F90`` variable in your makefile to be
an MPI compiler (e.g. ``mpifort``) and add ``-DMPI`` to the ``MPIFLAGGS``
variable:

.. code-block:: Makefile

    F90 = mpifort
    MPIFLAGS = -DMPI ... omitting other flags

However, you can also generate an MPI makefile by running the configuration
scripts and passing it MPI as the type.

.. important::

    When running ModEM with MPI you must use at least 2 tasks and
    a max of (2 x nTransmitters) + 1 tasks (except for the SP2 version compiled
    with -DFG).


Running ModEM: The basics
-------------------------

Running the compiled program with no command line arguments will print a brief
summary of usage

.. code-block:: text

    Output information to files, and progress report to screen (default).
    Usage: Mod3DMT -[job] [args]

    [READ_WRITE]
        -R rFile_Model rFile_Data [wFile_Model wFile_Data]
        Reads your input files and checks them for validity;
        optionally also writes them out
    [FORWARD]
        -F rFile_Model rFile_Data wFile_Data [wFile_EMsoln rFile_fwdCtrl]
        Calculates the predicted data and saves the EM solution
    [COMPUTE_J]
        -J rFile_Model rFile_Data wFile_Sens [rFile_fwdCtrl]
        Calculates and saves the full J(acobian)
    [MULT_BY_J]
        -M rFile_Model rFile_dModel rFile_Data wFile_Data [rFile_fwdCtrl]
        Multiplies a model by J to create a data vector
    [MULT_BY_J_T]
        -T rFile_Model rFile_Data wFile_dModel [rFile_fwdCtrl]
        Multiplies a data vector by J^T to create a model
    [INVERSE]
        -I NLCG rFile_Model rFile_Data [lambda eps]
        Here, lambda = the initial damping parameter for inversion
            eps = misfit tolerance for the forward solver
        OR
        -I NLCG rFile_Model rFile_Data [rFile_invCtrl rFile_fwdCtrl]
        Optionally, may also supply
            the model covariance configuration file [rFile_Cov]
            the starting model parameter perturbation [rFile_dModel]
        Runs an inverse search to yield an inverse model at every iteration
    [TEST_COV]
        -C rFile_Model wFile_Model [rFile_Cov]
        Applies the model covariance to produce a smooth model output
    [TEST_ADJ]
        -A J rFile_Model rFile_dModel rFile_Data [wFile_Model wFile_Data]
        Tests the equality d^T J m = m^T J^T d for any model and data.
        Optionally, outputs J m and J^T d.
    [TEST_SENS]
        -S rFile_Model rFile_dModel rFile_Data wFile_Data [wFile_Sens]
        Multiplies by the full Jacobian, row by row, to get d = J m.
        Compare to the output of [MULT_BY_J] to test [COMPUTE_J]

    Optional final argument -v [debug|full|regular|compact|result|none]
    indicates the desired level of output to screen and to files.


The ``-F [FORWARD]`` and ``-I [INVERSE]`` options will suffice for simple applications.
For some of the more advanced applications discussed below the ``-C [TEST_COV]`` option is
also needed. The ``-R [READ_WRITE]`` option can be useful for testing whether your input
files are in the correct format. The ``-v`` option controls the level of program output;
the default is “regular”; if you want to reduce outputs (“none” is *not recommended!*) you
can experiment with other output levels; “full ” and “debug” would be primarily useful for
developers. All of the other options (``-J``, ``-M``, etc.) are primarily useful for
development, testing, and some specialized applications that are beyond the scope of this
basic documentation.

Sub-directory ``examples/`` contains some sample input files with instructions for simple
cases. These should be sufficient to test your installation, and also will provide a
template for your own initial use. The ``doc/`` sub-directory contains some further
information on the ModEM system, beyond the summary of usage provided here.

All of the various options have a required set of arguments, which define input and output
file names (``rFile``, and ``wFile`` in the usage summary, respectively). There are also
optional arguments, mostly also specifying additional input or output file names. In all
cases the order of arguments must be exactly as indicated in the “Usage” block, so to give
an optional argument at the end of the list, all previous optional arguments must also be
provided. Thus, for example, with the ``-F`` option, if you want to specify the input
``rFile_fwdCtrl`` file, you must also specify the output ``wFile_EMsoln`` file, even if
you do not specifically want the full set of output EM fields, which this option will
result in. Before providing a more detailed discussion of the various options, we consider
file formats.

Using Docker to Compile and Run ModEM
-------------------------------------

.. note::

    For more information on Docker, please see: https://www.docker.com/resources/what-container/.

For convenience, a Dockerfile has been provided which will automatically build
ModEM once it is built. Once built, you can use the docker container to run
ModEM.  To run the ModEM docker container, first ensure you have `docker
installed <Docker-Getting-Started_>`_ and running (Either have the docker daemon
`docker running <daemon-running_>`_, or Docker `Desktop Running
<Docker-desktop_>`_.

The Dockerfile builds an Ubuntu instance with GFortran, MPI, and all the other
tools needed to compile and run ModEM and it will compile MF, SP, and SP2
versions of ModEM automatically.

Building
^^^^^^^^

.. warning::

    Before running ``docker build`` or ``docker run`` `the docker daemon <daemon-running_>`_
    or `Docker Desktop <Docker-desktop_>`_ must be running!

After you install Docker, and have either the docker deamon running (linux) or
docker desktop running (Mac/Windows) you can build the docker container using
the Dockerfile. Inside the ``ModEM`` directory (where the ``Dockerfile``
resides), run the following:

.. code-block:: bash

    $ cd ModEM
    $ docker build . -t modem:latest

After running the above command Docker will create an image that is taged with:
``modem:latest`` which can you run; it will have all necessary libraries and
compilers to compile and run ModEM.

As well, the Dockerfile will copy the ModEM code that it is in and place it into
`/home/modem/ModEM` in the Docker container and will build MF, SP and SP2
executables.

.. note::

    Building the ModEM Dockerfile the first time takes 2-4 minutes. After
    building once, Docker will cache the installed programs and later builds
    will take substantially less time.

Running ModEM in Docker
^^^^^^^^^^^^^^^^^^^^^^^^^

Now that we have built the docker container above, we can use ``docker run`` to
run the container and mount our local code:

.. code-block:: bash

    $ docker run --mount=type=bind,source=/absolute/path/to/ModEM,target=/home/modem/ModEM -it modem:latest

Normally, when running a docker container, the data created inside the container
will be deleted when the container terminates. To prevent this, we can use a
mount.

Mounting will allow your local ModEM code direcotry to be mounted into the
docker container. That way any change you make in one will be reflected on the
other.

After you terminate that container, and wish to start it again, you will need to
specify your run command with the same arguments as above to set the same mount.

The docker run command will run the image ``modem:latest`` and create a
container for us to use The ``-it`` argument will tell Docker to start an
interactive tty session and will drop you inside a bash shell inside the
running container.

Now, inside the docker container, we can move into ``~/ModEM/f90`` where we will
see three compiled ModEM executables:

* ``Mod3DMT_MF`` - matrix free version
* ``Mod3DMT_SP`` - SP version
* ``Mod3DMT_SP2`` - SP2 version

You can use these executables to run a ModEM example (which you will need to
mount into the Dockerfile yourself, see :ref:`docker-mounts`). Or, if you wish,
you an compile your own executable:

.. code-block:: bash

    $ cd ~/ModEM/f90
    $
    $ # Compile with MPI, SP2 solver:
    $ ./CONFIG/configure Makefile gfortran
    $ make
    $
    $ # or compile serial, SP2 solver
    $ ./CONFIG/configure -m serial Makefile gfortran
    $ make clean
    $ make
    $
    $ # Run ModEM serially or with MPI:
    $ mpiexec -n 2 ./Mod3DMT ...
    $ # Or if you compiled serial:
    $ ./Mod3DMT ...

You can terminate this docker session by logging out by typing ``exit`` or you
can use ``Ctrl-D``.

.. important::

    Again, if you want to re-run docker and reuse the files/folders/data you
    created before you will need to specify the volume that you used before

    .. code-block:: bash

        $ docker run --mount=type=bind,source=/absolute/path/to/ModEM,target=/home/modem/ModEM -it modem:latest

.. warning::

    f you want to save data you will need to specify additional volumes or
    mounts. See :ref:`docker-mounts`.

    For more information on persistent data inside a container see:

    https://docs.docker.com/get-started/docker-concepts/running-containers/persisting-container-data/

.. _docker-mounts:

Specifying mounts points
^^^^^^^^^^^^^^^^^^^^^^^^^

.. note:: 

    For more information on Docker mounts see:
    https://docs.docker.com/engine/storage/bind-mounts/.

You may want to add additional mounts in your docker container. For example, you
might want to include the ModEM-Examples direcotry:

.. code-block:: bash

    $ docker run --mount=type=bind,source=/absolute/path/to/ModEM,target=/home/modem/ModEM \
                 --mount=type=bind,source=/absolute/path/to/ModEM-Examples,target=/home/modem/ModEM-Examples \
                 -it modem

.. _Docker-desktop: https://www.docker.com/blog/getting-started-with-docker-desktop/
.. _daemon-running: https://docs.docker.com/engine/daemon/start/
.. _docker-Getting-Started: https://www.docker.com/get-started/
.. _dockerfile: https://github.com/magnetotellurics/ModEM/blob/main/Dockerfile

.. important::

    Again, each time you run your modem container, you will need to specify
    volume and mount points.
