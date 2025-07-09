Global3D / Earth ModEM
======================

This version of ModEM is a global version that uses a forward solver that was
created by Uyeshima & Schultz (2000). You can compile it in a similar way to
normal ModEM:

```
$ ./CONFIG/Configure.EARTH.OSU.GFortran Makefile Release
```

The configuration files have some difficulty creating MPI Makefiles, except for
`./CONFIG/Configure.EARTH.OSU.Intel`. You can use the `Makefile.gnu.mpi` to
compmile with MPI:

```
$ make -f Makefile.gnu.mpi
```

An example of a Global run can be found in [ModEM-Examples][ModEM-Examples]

[ModEM-Examples]: https://github.com/MiCurry/ModEM-Examples


If you use this branch in new research or derivative work, please include the 
following citation:

M. Uyeshima, A. Schultz, Geoelectromagnetic induction in a heterogeneous
sphere:a new three-dimensional forward solver using a conservative
staggered-grid finite difference method, Geophysical Journal International,
Volume 140, Issue 3, March 2000, Pages 636–650,
https://doi.org/10.1046/j.1365-246X.2000.00051.x
