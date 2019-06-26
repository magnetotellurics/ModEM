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
rm -rf outputs/test_forward/
#
# CREATE TEST OUTPUT FOLDER
mkdir -p test_forward
#
# ENTER TEST OUTPUT FOLDER
cd test_forward/
#
#
echo "#### START FORWARD MPI TEST WITH $ncores CORES AT $now ####\n" >> std_out.txt
#
#
echo "#### COMMAND LINE: 'mpirun -n $ncores ../$EXEC -R ../inputs/rFile_Model.sw ../inputs/rFile_Data.dat -v full'\n" >> std_out.txt
#
#
mpirun -n $ncores ../$EXEC -F ../inputs/rFile_Model.sw ../inputs/rFile_Data.dat ../inputs/wFile_Data.dat -v full &>> std_out.txt
#
#
echo "### FINISH FORWARD MPI TEST ###\n" >> std_out.txt
#
#
cd ..
#
#
mv test_forward/ outputs/
#
#
# END OF SCRIPT

