#!/bin/bash
#
# ARGUMENTS: 1 - Mod3DMT EXECUTABLE, 2 - MODEL FILE, 3 - DATA FILE
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
#
# STRING NOW
NOW=$(date "+%Y/%m/%d - %H:%M:%S")
#
# GET ENVIROMENT NUMBER OF CORES 
NCORES=$(nproc)
#
# CREATE TEST OUTPUT FOLDER
mkdir -p ${EXEC_NAME}
#
# ENTER TEST OUTPUT FOLDER
cd ${EXEC_NAME}/
#
#
echo "#### START MULT_BY_J_T MPI TEST WITH $NCORES CORES AT $NOW ####" | tee std_out.txt
#
#
echo "#### COMMAND LINE: [mpirun -n $NCORES ../$EXEC -T ../$MODEL ../$DATA wFile_dModel -v full]" | tee -a std_out.txt
#
#
mpirun -n $NCORES ../$EXEC -T ../$MODEL ../$DATA wFile_dModel -v full | tee -a std_out.txt
#
# CATCH RESULT
result=$?
#
# TEST RESULT
if [ "$result" -ne "0" ]; then
	#
	#
	echo "TEST MULT_BY_J_T FAIL: $result" | tee -a std_out.txt
	#
	#
	cd ..
	#
	#
	exit $result
fi
#
#
echo "#### FINISH MULT_BY_J_T MPI TEST ####" | tee -a std_out.txt
#
#
cd ..
#
#
mv ${EXEC_NAME}/ outputs/temp/test_mult_by_j_t
#
#
exit 0
#
# END OF SCRIPT

