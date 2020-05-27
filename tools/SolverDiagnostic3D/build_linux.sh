#!/bin/bash
#
mkdir -p bin
#
cd bin/
#
#g++ -c ../src/*.h
#
g++ ../src/*.cpp -O3 -o SolverDiagnostic3D
#
cd ..
#

