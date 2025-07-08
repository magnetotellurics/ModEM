# Modular Electromagnetic Inversion Software (ModEM) - CSEM Version

This version of ModEM contains code to run Controlled Source Electro-Magnetic
forward and inversions and was developed by Brazil's Observatório Nacional
(ON).

We've included the CSEM branch as an experimental version of ModEM. There are a
few CSEM examples in the [ModEM-Examples repository][ModEM-Examples] that run
with EM1D (See note). This branch should serve as a good starting place
for future CSEM work in ModEM.

> **NOTE:** EM1D is currently not included in ModEM-Model due to reasons
> described below in [CSEM Forward Programs](#csem-forward-programs).

## Related Repositories

* [ModEM-Tools][ModEM-Tools] - A collection of MatLab and Python tools
to manipulate ModEM input and output files.
* [ModEM-Examples][ModEM-Examples] - A collection of 2D and 3D
examples

[ModEM-Tools]: https://github.com/MiCurry/ModEM-Tools
[ModEM-Examples]: https://github.com/MiCurry/ModEM-Examples

## CSEM Forward Programs

This branch can run with two CSEM codes, EM1D, created by Rita Streich, and
[Dipole1D][dipole1D] crated by Scripps' Marine EM Laboratory.

Currently, EM1D runs the best of ModEM and the examples found in 
[ModEM-Examples][ModEM-Examples]; however, the EM1D code was removed from this
repository as it contained propriety code that we cannot include. Work is being
done to remove these parts of EM1D and we hope to include an updated EM1D to 
this repository.


Dipole1D does not run as well with the examples in
[ModEM-Examples][ModEM-Examples]. It still might server as a good starting place
for future work with CSEM and perhaps better examples can be created.

We are not able to include Dipole1D in this repository due to it's license; however,
it is still free and open source. To use Dipole1D, you can either:

1. Use the `src/CONFIG/Configue` script to automatically download and extract
Dipole1D - or -
2. Manually download [OCCAM1DCSEM][dipole1d] and extract it and place the .f90
files inside the `source/` directory inside `src/3D_MT/CSEM_Module/Dipole1D/`
directory.

[dipole1d]: https://marineemlab.ucsd.edu/Projects/Occam/1DCSEM/

## Compiling CSEM Branch

### Configuring

Unlike ModEM, the CSEM branch contains a single CSEM configuration file. We can
use it to create a CSEM capable version of ModEM. There are two CSEM forward
solvers available with the ModEM CSEM branch [Dipole1D][dipole1d] and EM1D.

EM1D is included within the CSEM branch, but Dipole1D is not due to Dipole1D's
license. Thus, Dipole1D will need to be downloaded before it is run,
thankfully the Configuration script can automatically download and extract
Dipole1D into the correct location. See [Dipole1D
Configuration](#dipole1d-configuration)

[dipole1d]: https://marineemlab.ucsd.edu/Projects/Occam/1DCSEM/

#### Configuration Script options

Running the Configuration script without arguments will produce the usage
message:

```bash
$ ./CONFIG/Configure
Usage: ./CONFIG/Configure with the following options:
Compiler: Choose from supported compilers: [ gfortran | ifort ]
Makefile: Provide a name for your output Makefile name.
[Debug or Release]: Choose whether you want to compile the Debug or Release version.
[MPI or Serial]:  Choose whether you want to compile the parallel (MPI) or
 serial version.
[MF or SP or SP2]:  Choose between the Matrix Free (MF), or the Modified System
    of Eqs 1 (SP), or the Modified System of Eqs 2 (SP2) of the code.
[MT or MT+CSEM]:  Compile MT or MT+CSEM. In Case of MT+CSEM, choose in the
    following option whether Dipole1D or EM1D or both will be used to get for the
    secondary field formulation.
[Dipole1D or EM1D or Dipole1D+EM1D]:  (Optional) - Choose whether you have
    Dipole1D, or EM1D or both codes in the source files folder '/3D_MT/CSEM_module'
Optional: Enviornment variables: 'FC' 'FFLAGS' 'CPPFLAGS' 'LDFLAGS' 'LDLIBS'
are respected
```

#### Dipole1D Configuration

To configure a ModEM executable that uses Dipole1D, you can use the following
configuration:

```bash
$ cd ModEM-Model/src/
$ ./CONFIG/Configure gfortran Makefile Release Serial SP2 MT+CSEM Dipole1D
Dipole1D is not currently in ./3D_MT/CSEM_module/Dipole1D/
Would you like to have this script automatically download it now? [Yes/No]:
```

If you don't have Dipole1D downloaded and extracted in the
`src/3D_MT/CSEM_module/Dipole1D` directory you will get the above message.
Passing 'Yes/yes/Y/y' will automatically download and extract Dipole1D into the
correct location.

If there are problems with the download script, or you would like to download
Dipole1D by yourself, you can download Dipole1D yourself and extract the source
code into: `src/3D_MT/CSEM_module/Dipole1D`. (i.e. ensure that `Dipole1D.f90`
and other Dipole1D.f90 files are in `src/3D_MT/CSEM_module/Dipole1D`).

After downloading CSEM (Via the script), you'll get the following message:

```bash
# Using compile cmd gfortran from cmd line
# Using compiler optimization options -cpp -DCSEM_DIPOLE1D -O3
# -ffree-line-length-none -dI -fallow-argument-mismatch from cmd line
# Using compiler MPI flags -cpp -DCSEM_DIPOLE1D from cmd line
# Using ./objs/3D_MT/csemBuild for object file output directory
# Using Link options -lblas -llapack from cmd line
# Using Library path  from cmd line
# Using search path from cmd line:
 .:MPI:INV:SENS:UTILS:FIELDS:FIELDS/FiniteDiff3D:3D_MT:3D_MT/DICT:3D_MT/ioMod:3D_MT/modelParam:3D_MT/modelParam/modelCov:3D_MT/modelParam/modelParamIO:3D_MT/CSEM_module:3D_MT/FWD_SP2:3D_MT/SP_Topology:3D_MT/FWD:3D_MT/FWD/Mod2d:3D_MT/CSEM_module:3D_MT/CSEM_module/Dipole1D
 Couldn't find source file for module EM1D
```

The last warning message can be ignored because we did not request EM1D.

#### EM1D Configuration

> **NOTE:** EM1D is currently not included in ModEM-Model due to reasons
> described above in [CSEM Forward](#csem-forward-programs).

In order to run the EM1D version of ModEM, you will need to include the location
of the [FFTW library][fftw] installation. To do so you can pass the installation
in the `FFTW` environment variable:

```bash
export FFTW=/path/to/FFTW/
```

(The FFTW should point to the folder that contains both the `lib/` and the
`include/` directories.)

For more information on installing the FFTW library see: [FFTW
Installation](#fftw-install).

```
$ export FFTW=/path/to/FFTW_INSTALL/
$ ./CONFIG/Configure gfortran Makefile Release Serial SP2 MT+CSEM EM1D
```

[fftw]: https://www.fftw.org/

#### Dipole1D+EM1D Configuration

To run with both Dipole1D+EM1D, you will similarly you will need to set the
FFTW variable in order to run the EM1D:

```
$ export FFTW=/path/to/FFTW_INSTALL/
$ ./CONFIG/Configure gfortran Makefile Release Serial SP2 MT+CSEM Dipole1D+EM1D
```

#### Cray system?

If you're using a Cray System, you can use the FC variable during your
configuration script to change the name of the compiler to use for your
makefile:

```
$ export FC=ftn
$ ./CONFIG/Configure gfortran Makefile Release Serial SP2 MT+CSEM Dipole1D
```

#### Making 

Once our Makefile is created we can rename it to Makefile, if we named it
something other than Makefile and then call make:

```
make
```
