#!/bin/bash
#
chmod 777 Build/*
#
./Build/Configure.modem.ubuntu.gfortran MakefileMPI MPI ModEM.f90
make -f MakefileMPI clean
make -f MakefileMPI
mv ModEM ModEM_MPI
#