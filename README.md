Modular Electromagnetic Inversion Software (ModEM)
==================================================

**ModEM: A modular system for inversion of electromagnetic geophysical data.**  
**Authors**: Gary Egbert, Anna Kelbert, Naser Megbel & Hao Dong.

> **NOTE:** This repository has been converted from the ModEM's OSU CEOAS Subversion
> (SVN) repository. SVN revisions have been preserved and converted into Git
> commits. Some main branches have been renamed:
>
> | SVN Branch Name | GitHub Branch Name |
> | --------------- | -------------------|
> | trunk           | [trunk][trunk-branch] |
> | stable-         | [main][main-branch] |
> | stable          | [classic][classic-branch] |
>
> Furthermore, the `matlab` and `examples` directory have been moved into the
> [ModEM-Tools][ModEM-Tools] and [ModEM-Examples][ModEM-Examples] repositories,
> respectively.

[trunk-branch]: https://github.com/magnetotellurics/ModEM/tree/trunk
[main-branch]: https://github.com/magnetotellurics/ModEM/tree/main
[classic-branch]: https://github.com/magnetotellurics/ModEM/tree/classic

# Contents

* [Obtaining The Software](#obtaining-the-software)
* [Building ModEM](#building-modem)
    * [Dependencies](#dependencies)
    * [Creating Makefiles From Configuration files](#creating-makefiles-from-configuration-files)
    * [Compiling](#compiling)
    * [Compiling with MPI](#compiling-with-mpi)
* [Basic ModEM Usage](#basic-modem-usage)
    * [Forward Modeling](#forward-modeling)
    * [Inversion Modeling](#inversion-modeling)
* [Building and Running in Docker](#building-and-running-in-docker)
* [More Information And Tools](#more-information-and-tools)
    * [Related Repositories](#related-repositories)
    * [ModEM-ON and ModEM-OO](#modem-on-and-modem-oo)
* [Citations](#citations)

# Obtaining the Software 

You can download this source code by cloning or downloading this repository. To
clone with git run `git clone`:

``` bash
$ git clone https://github.com/magnetotellurics/ModEM.git
Cloning into 'ModEM'...
remote: Enumerating objects: 18608, done.
remote: Counting objects: 100% (199/199), done.
remote: Compressing objects: 100% (54/54), done.
remote: Total 18608 (delta 162), reused 152 (delta 145), pack-reused 18409 (from 1)
Receiving objects: 100% (18608/18608), 73.43 MiB | 13.34 MiB/s, done.
Resolving deltas: 100% (14565/14565), done.
```

For more information on Git, GitHub and cloning, please see:
https://docs.github.com/en/get-started.

You can also download specific versions and commits of ModEM by following [these
instructions][GitHub-Download-Tutorial] provided by GitHub. Please note though
that these downloads do not contain any git repository history or information.

[GitHub-Download-Tutorial]:https://docs.github.com/en/repositories/working-with-files/using-files/downloading-source-code-archives

# Building ModEM

## Dependencies

ModEM depends on the [LAPACK][lapack] and [BLAS][blas] libraries. Often, these
are already included on most systems, but if not you will need to install them
yourself and ensure they are properly linked in the `LIBS` and `LIBS_PATH`
variables of the Makefile.

[lapack]:https://www.netlib.org/lapack/
[blas]:https://www.netlib.org/blas/

## Creating Makefiles from Configuration files

The current build system for ModEM uses the `fmkmf.pl` pearl script and the
configuration scripts that are in `f90/CONFIG`. Although ModEM uses Make and
Makefiles, makefiles are meant to be created by these configuration scripts.

To create Makefiles, run the configuration scripts that match your system,
desired compiler and desired ModEM version. These configuration scripts
call the `fmkmf.pl` script:

``` bash
$ cd f90/
$ ./CONFIG/Configure.3D_MT.MAC.GFortran makefile.gnu.mpi MPI
```

All configuration scripts take the same arguments:

```bash
$ /CONFIG/Configfile <desired-makefile-name> <type>
```

Where `type` is either: `MPI, release, debug`:

* `MPI` - Generates a makefile that will compile an MPI version of ModEM
* `release` - Generates a makefile that will compile a serial version of ModEM
* `debug` -  Generates a makefile that will compile a serial version of ModEM
 with `-O0` (no optimizations) and debug symbols (`-g`, bounds checking ,etc)

Of course, any generated makefile can be altered as you see fit. For instance,
it is often helpful to add the `-g` to copmile with debug symbols when working
with an MPI makefile.

## Compiling

Once you have generated your own makefile, or decided to use one of the
defaults, you can compile ModEM by renaming that makefile to `Makefile` and
then run `make`:

``` bash
$ cp Makefile.3D.MF.gnu Makefile
$ make clean # Not necessary from a fresh clone
$ make
```

### Compiling with MPI

If you have MPI installed, and already have a serial makefile generated you can
easily compile the MPI version of ModEM by altering the makefile. This is
sometimes easier then creating a new makefile from the configuration scripts.

To compile an MPI version, edit the compiler in your makefile to be an MPI
compiler (e.g.  `mpifort`) and add `-DMPI` to the `MPIFLAGGS` variable:

``` Makefile
F90 = mpifort
MPIFLAGS = -DMPI ... omitting other flags
```

However, you can also generate an MPI makefile by running the configuration scripts
and passing it `MPI` as the type.

> **IMPORTANT NOTE:** When running ModEM with MPI you must use at least 2 tasks and
> a max of `(2 x nTransmitters) + 1` tasks (except for the SP2 version
> compiled with `-DFG`).

# Basic ModEM Usage

> **NOTE:** For a more detailed information on ModEM usage please see the [ModEM
> User's Guide][Users-Guide].

While the [User's Guide][Users-Guide] is the best resource for information on
running ModEM, information can also be found by running the `Mod2DMT` or
`Mod3DMT` executables with no arguments. Furthurmore, specifying the job flag
with no other arguments will produce a detailed usage description for that job
type.  For example: `./Mod3DMT -F` will produce a detailed description of
options for forward modeling.

The examples below will use data and model files from the
[BLOCK2][BLOCK2-Example] Magnetotelluric examples found within the
[ModEM-Examples][ModEM-Examples] repo. You will need to clone or download this
repository to obtain the data files and examples.

[BLOCK2-Example]: https://github.com/magnetotellurics/ModEM-Examples/tree/main/Magnetotelluric/3D_MT/BLOCK2
[Users-Guide]: https://github.com/magnetotellurics/ModEM/blob/main/doc/userguide/ModEM_UserGuide.pdf

## Forward Modeling

Once we have ModEM compiled, we can link or copy the executable into our BLOCK2
example:

``` bash
$ cd ModEM-Examples/Magnetotelluric/3D_MT/BLOCK2
$ ln -s ~/ModEM/f90/Mod3DMT . # Link or copy if you prefer.
```

Then we can run ModEM with `-F` to run the forward model:

```bash
$ ./Mod3DMT -F m1.ws Templdate_de.dat fwd.BLOCK2.dat esoln.BLOCK2.dat
```

This will calculate the predicted data from the model `m1.ws` for the periods,
station locations and data types in `Templdate_de.dat` and will write out of the
results in `fwd.BLOCK2.dat`. The last argument `esoln.BLOCK2.dat` is an optional
argument and if present will write out the full electro-magnetic solution in 
the specified file.

If you compiled ModEM with MPI, you can run the above example with MPI:

```bash
$ mpiexec -n 2 ./Mod3DMT -F m1.ws Templdate_de.dat fwd.BLOCK2.dat esoln.BLOCK2.dat
```

When running ModEM with MPI, you **must** use at least **two** MPI tasks. For
all versions of ModEM, except SP2 compiled with `-DFG`, the max number of MPI
tasks you can use is `(2 x nTransmitters) + 1`. Where the number of transmitters
is the number of periods/frequencies. Thus you can run with 9 tasks: 4
transmitters multipled by 2 polarizations + 1 main task:

```bash
$ mpiexec -n 9 ./Mod3DMT -F m1.ws Templdate_de.dat fwd.BLOCK2.dat esoln.BLOCK2.dat
```

Of course, this will work only if you have 9 cores available on your system.

> **IMPORTANT NOTE:** Running ModEM with MPI you must use at least 2 tasks and
> a max of `(2 x nTransmitters) + 1` tasks (except for the SP2 version
> compiled with `-DFG`).

> **NOTE:** Please see the [User's Guide][Users-Guide] for additional, optional arguments for
> forward modeling.

## Inversion Modeling

We can run the inverse search to create a model from data with the following arguments 
for our BLOCK2 example:

``` bash
$ ./Mod3DMT -I NLCG pr.ws de.dat block2.inv
```

This will search for a model, starting from `pr.ws` and that matches the
dimensions and size of `pr.ws`, from the data in `de.dat` using the non-linear
conjugate gradient (NLCG).

Options after `de.dat` are optional, and can be either of the following:

```
-I NLCG rFile_Model rFile_Data [lambda eps]
or
-I NLCG rFile_Model rFile_Data [rFile_invCtrl rFile_fwdCtrl]
or
I NLCG rFile_Model rFile_Data [InvCtrl FwdCtrl CovCtrl StartModel]
```

| Options | Description | 
| --------| ------------| 
| lambda | initial damping parameter for inversion |
| eps | misfit tolerance for the forward solver |
| rFile_invCtrl | Inversion control file |
| rfile_fwdCtrl | Forward control file | 
| CovCtrl |  Covariance file |
| StartModel | Starting model file |

For more information on these arguments, and formats of the control files, and
covariance files please see the [ModEM User's Guide][Users-Guide].

Similar to the forward modeling, you can run the inverse modeling with MPI. The
same rules as apply as to the number of transmitters and tasks: Must use at least two
tasks and no more than (2 x nTransmitters) + 1 tasks (except for the SP2 version
compiled with `-DFG`):

``` bash
$ mpiexec -n 9 ./Mod3DMT -I NLCG pr.ws de.dat block2.inv
```

# Building and Running In Docker

A [Dockerfile][Dockerfile] has been provided that will create an Ubuntu docker
container where you can compile and run ModEM. To run the ModEM docker
container, first ensure you have [docker installed][Docker-Getting-Started] and
running (Either have the docker [deamon running][Damon-running], or Docker
[Desktop running][Docker-desktop].

## Building

After you install Docker, you will need to build the docker container. Inside
the `ModEM` directory (where the `Dockerfile` resides), run the following:

```bash
$ docker build . -t modem:latest --build-args ncpus=2
```

You can specify different numbers of arguments in the `ncpus` build argument if
you wish. It will speed up the time when Docker compiles MPICH.

After running the above command Docker will create an image named
`modem:latest` that you can then use to run; it will have all necessary
libraries and compilers to compile and run ModEM.

## Running and building ModEM in Docker

To run the container, we need to run `docker run`; however, we will also need
to pass in a `--mount` option to mount our source code into our docker
container:

```bash
$ docker run --mount=type=bind,source=/abosolute/path/to/ModEM,target=/root/ModEM -it modem
```

This will run the container and the `-it` argument will tell Docker to start an
interactive tty session and will drop you inside the running container. You
will find that the ModEM source code that you specified in the `source` mount
argument will be present. Any changes you make in this folder, in either the
container or the host machine will be present on the other machine.

Once you're inside the container, `cd` into `ModEM/f90`. There, you can run the
`./CONFIG/Configure.3D_MT.Docker.GFortran` and `make` as described in other
sections above.

## Specifying additional mounts

It might be helpful to specify additional mount points, for instance it might
be helpful to mount a work directory, such as the ModEM-Examples repository:

`
```bash
$ docker run --mount=type=bind,source=/abosolute/path/to/ModEM,target=/root/ModEM \
             --mount=type=bind,soucr=/home/users/uname/ModEM-Examples,target=/root/ModEM-Examples \
             -it modem
```

[Docker-desktop]: https://www.docker.com/blog/getting-started-with-docker-desktop/
[Damon-running]: https://docs.docker.com/engine/daemon/start/
[Docker-Getting-Started]: https://www.docker.com/get-started/
[Dockerfile]:./Dockerfile

# More Information and Tools

## Related Repositories

* [ModEM-Tools][ModEM-Tools] - A collection of MatLab and Python tools
to manipulate ModEM input and output files.
* [ModEM-Examples][ModEM-Examples] - A collection of 2D and 3D
examples

[ModEM-Tools]: https://github.com/magnetotellurics/ModEM-Tools
[ModEM-Examples]: https://github.com/magnetotellurics/ModEM-Examples

## ModEM-ON and ModEM-OO

The versions of ModEM worked on by Brazil's Observatório Nacional (ON) are
included as branches in this repository under the following names:

| ON Repo Name | GitHub Branch name |
| ------------ | ------------------ |
| ModEM-OO Repository | [object-oriented][Object-Oriented-Branch] |
| ModEM-ON Repository | [CSEM][CSEM-Branch] |

> <font color='red'>**WARNING**:</font> Both the object-oriented and CSEM
> branch have unrelated git-history than other branches. Thus, these branches
> should not be merged into other branches.

[Object-Oriented-Branch]: https://github.com/magnetotellurics/ModEM/tree/object-oriented
[CSEM-Branch]: https://github.com/magnetotellurics/ModEM/tree/CSEM

# Citations

If you use ModEM in your work, please cite the following two sources:

Anna Kelbert, Naser Meqbel, Gary D. Egbert, Kush Tandon, ModEM: A modular system
for inversion of electromagnetic geophysical data, Computers & Geosciences,
Volume 66, 2014, Pages 40-53, ISSN 0098-3004,
https://doi.org/10.1016/j.cageo.2014.01.010.

Egbert, Gary & Kelbert, Anna. (2012). Computational recipes for electromagnetic
inverse problems. Geophysical Journal International. 189. 251-267.
10.1111/j.1365-246X.2011.05347.x. 

If you use the SP versions (divergence-free forward method) in your work, please 
cite this source  as well: 

Hao Dong,  Gary D. Egbert. Divergence-free solutions to electromagnetic forward 
and adjoint problems: a regularization approach, Geophysical Journal International,
Volume 216(2), 2019, Pages 906–918, https://doi.org/10.1093/gji/ggy462

If you use the CUDA/HIP interfaces of ModEM in your work, please also
cite the following source:

Hao Dong, Kai Sun, Gary Egbert, Anna Kelbert, Naser Meqbel, Hybrid CPU-GPU
solution to regularized divergence-free curl-curl equations for electromagnetic
inversion problems, Computers & Geosciences, Volume 184, 2024, 105518, ISSN
0098-3004, https://doi.org/10.1016/j.cageo.2024.105518.
