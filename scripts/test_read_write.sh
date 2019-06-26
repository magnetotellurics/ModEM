#!/bin/bash
#
# ARGUMENTS: 1 - Mod3DMT EXECUTABLE
EXEC=$1
#
# STRING NOW
now=$(date "+%Y/%m/%d - %H:%M:%S")
#
# GET ENVIROMENT NUMBER OF CORES 
ncores=$(nproc)
#
# CREATE TEST OUTPUT FOLDER
mkdir -p test_read_write
#
# ENTER TEST OUTPUT FOLDER
cd test_read_write/
#
#
echo "#### START READ_WRITE MPI TEST WITH $ncores CORES AT $now ####" >> std_out.txt
#
#
mpirun -n $ncores ../$EXEC ../inputs/rFile_Model.sw ../inputs/rFile_Data.dat &>> std_out.txt
#
#
echo "### FINISH READ_WRITE MPI TEST ###" >> std_out.txt
#
#
cd ..
#
#
mv test_read_write/ outputs/
#
#
# END OF SCRIPT

