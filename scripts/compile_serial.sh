#!/bin/bash
#
./Configure.modem.ubuntu.gfortran MakefileSerial SERIAL ModEM.f90
make -f MakefileSerial clean
make -f MakefileSerial
mv ModEM ModEM_SERIAL
#