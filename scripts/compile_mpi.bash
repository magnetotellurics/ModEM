#!/bin/bash
#
./Configure.modem.ubuntu.gfortran MakefileMPI MPI ModEM.f90
make -f MakefileMPI clean
make -f MakefileMPI
mv ModEM ModEM_MPI
#