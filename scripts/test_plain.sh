#!/bin/bash
#
# ARGUMENTS: 1 - Mod3DMT EXECUTABLE, 2 - NUMBER_OF CORES
EXEC=$1
EXEC_NAME=$EXEC
emptyspace=""
EXEC_NAME=${EXEC_NAME/.exe/$emptyspace}
EXEC_NAME=${EXEC_NAME/.sh/$emptyspace}
EXEC_NAME=${EXEC_NAME/.txt/$emptyspace}
EXEC_NAME=${EXEC_NAME/*\//$emptyspace}
#
ncores=$2
#
# STRING NOW
now=$(date "+%Y/%m/%d - %H:%M:%S")
#
# GET ENVIROMENT NUMBER OF CORES 
ncores=$(nproc)
#
# REMOVE TEST OUTPUT FOLDER FROM MAIN OUTPUT FOLDER
rm -rf outputs/test_plain_${EXEC_NAME}/
#
# CREATE TEST OUTPUT FOLDER
mkdir -p test_plain_${EXEC_NAME}
#
# ENTER TEST OUTPUT FOLDER
cd test_plain_${EXEC_NAME}/
#
#
echo "#### START MODEM MPI TEST WITH $ncores CORES AT $now ####" | tee -a std_out.txt
#
#
echo "#### COMMAND LINE: [mpirun -n $ncores ../$EXEC -v full]" | tee -a std_out.txt
#
#
mpirun -n $ncores ../$EXEC -v full | tee -a std_out.txt
#
# CATCH RESULT
result=$?
#
# TEST RESULT
if [ "$result" -ne "0" ]; then
	#
	#
	echo "TEST PLAIN FAIL: $result" | tee -a std_out.txt
	#
	#
	cd ..
	#
	#
	exit $result
fi
#
#
echo "#### FINISH PAIN MPI TEST ####" | tee -a std_out.txt
#
#
cd ..
#
#
mv test_plain_${EXEC_NAME}/ outputs/temp/
#
#
exit 0
#
# END OF SCRIPT

