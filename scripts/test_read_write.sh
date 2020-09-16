#!/bin/bash
#
# ARGUMENTS: 1 - Mod3DMT EXECUTABLE, 2 - MODEL FILE, 3 - DATA FILE, 4 - NUMBER_OF CORES
EXEC=$1
EXEC_NAME=$EXEC
emptyspace=""
EXEC_NAME=${EXEC_NAME/.exe/$emptyspace}
EXEC_NAME=${EXEC_NAME/.sh/$emptyspace}
EXEC_NAME=${EXEC_NAME/.txt/$emptyspace}
EXEC_NAME=${EXEC_NAME/*\//$emptyspace}
#
MODEL=$2
DATA=$3
ncores=$4
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
echo "#### START READ_WRITE MPI TEST WITH $ncores CORES AT $now ####" | tee -a ${EXEC_NAME)_std_out.txt
#
#
echo "#### COMMAND LINE: [mpirun -n $ncores ../$EXEC -R ../$MODEL ../$DATA wFile_Model.sw wFile_Data.dat -v full]" | tee -a ${EXEC_NAME)_std_out.txt
#
#
mpirun -n $ncores ../$EXEC -R ../$MODEL ../$DATA wFile_Model.sw wFile_Data.dat -v full | tee -a ${EXEC_NAME)_std_out.txt
#
# CATCH RESULT
result=$?
#
# TEST RESULT
if [ "$result" -ne "0" ]; then
	#
	#
	echo "TEST READ_WRITE FAIL: $result" | tee -a ${EXEC_NAME)_std_out.txt
	#
	#
	cd ..
	#
	#
	exit $result
fi
#
#
echo "#### FINISH READ_WRITE MPI TEST ####" | tee -a ${EXEC_NAME)_std_out.txt
#
#
cd ..
#
#
mv test_read_write/  outputs/temp/
#
#
exit 0
#
# END OF SCRIPT

