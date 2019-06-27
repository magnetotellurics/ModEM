#!/bin/sh
EXEC=$1
echo $1
pwd
ncores=4
mpirun -n $ncores $EXEC -F rFile_Model.ws rFile_Data.dat wFile_Data.dat
exit 0

