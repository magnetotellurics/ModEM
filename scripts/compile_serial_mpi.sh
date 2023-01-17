#!/bin/bash
#
./Configure.modem.ubuntu.gfortran MakefileSerial SERIAL ModEM.f90
make -f MakefileSerial clean
make -f MakefileSerial
mv ModEM ModEM_SERIAL
#
./Configure.modem.ubuntu.gfortran MakefileMPI MPI ModEM.f90
make -f MakefileMPI clean
make -f MakefileMPI
mv ModEM ModEM_MPI
#