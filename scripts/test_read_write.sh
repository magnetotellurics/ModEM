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
# REMOVE TEST OUTPUT FOLDER FROM MAIN OUTPUT FOLDER
rm -rf outputs/test_read_write/
#
# CREATE TEST OUTPUT FOLDER
mkdir -p test_read_write
#
# ENTER TEST OUTPUT FOLDER
cd test_read_write/
#
#
echo "#### START READ_WRITE MPI TEST WITH $ncores CORES AT $now ####\n" >> std_out.txt
#
#
echo "#### COMMAND LINE: 'mpirun -n $ncores ../$EXEC -R ../inputs/rFile_Model.sw ../inputs/rFile_Data.dat -v full'\n" >> std_out.txt
#
#
mpirun -n $ncores ../$EXEC -R ../inputs/rFile_Model.sw ../inputs/rFile_Data.dat wFile_Model.sw /wFile_Data.dat -v full &>> std_out.txt
#
#
echo "### FINISH READ_WRITE MPI TEST ###\n" >> std_out.txt
#
#
cd ..
#
#
mv test_read_write/ outputs/
#
#
# END OF SCRIPT

