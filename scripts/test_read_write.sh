#!/bin/bash
#
# ARGUMENTS: 1 - Mod3DMT EXECUTABLE, 2 - MODEL FILE, 3 - DATA FILE
EXEC=$1
MODEL=$2
DATA=$3
#
# STRING NOW
now=$(date "+%Y/%m/%d - %H:%M:%S")
#
# GET ENVIROMENT NUMBER OF CORES 
ncores=$(nproc)
#
#
echo "#### START READ_WRITE MPI TEST WITH $ncores CORES AT $now ####\n" >> std_out.txt
#
#
echo "#### COMMAND LINE: [mpirun -n $ncores $EXEC -R $MODEL $DATA wFile_Model.sw wFile_Data.dat -v full]" >> std_out.txt
#
#
mpirun -n $ncores $EXEC -R $MODEL $DATA wFile_Model.sw wFile_Data.dat -v full &>> std_out.txt
#
# CATCH RESULT
result=$?
#
# TEST RESULT
if [ "$result" -ne "0" ]; then
	#
	#
	echo "TEST READ_WRITE FAIL: $result"
	#
	#
	exit $result
fi
#
#
echo "### FINISH READ_WRITE MPI TEST ###\n" >> std_out.txt
#
# REMOVE OLD TEST FOLDER RESULT FROM outputs/
rm -rf outputs/test_read_write
#
# CREATE NEW TEST FOLDER RESULT ON outputs/
mkdir outputs/test_read_write
#
# MOVE ALL FILES TO NEW TEST FOLDER RESULT ON outputs/
mv * outputs/test_read_write
#
#
# END OF SCRIPT

