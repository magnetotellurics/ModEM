#!/bin/sh
EXEC=$1
echo $1
pwd
ncores=4
mpirun -n $ncores ../$EXEC
exit 0

