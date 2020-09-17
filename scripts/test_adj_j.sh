#!/bin/bash
#
# ARGUMENTS: 1 - Mod3DMT EXECUTABLE, 2 - MODEL FILE, 3 - dMODEL FILE, 4 - DATA FILE
EXEC=$1
EXEC_NAME=$EXEC
emptyspace=""
EXEC_NAME=${EXEC_NAME/.exe/$emptyspace}
EXEC_NAME=${EXEC_NAME/.sh/$emptyspace}
EXEC_NAME=${EXEC_NAME/.txt/$emptyspace}
EXEC_NAME=${EXEC_NAME/*\//$emptyspace}
#
MODEL=$2
dMODEL=$3
DATA=$4
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
echo "#### START ADJ J MPI TEST WITH $NCORES CORES AT $NOW ####" | tee std_out.txt
#
#
echo "#### COMMAND LINE: [mpirun -n $NCORES ../$EXEC -A J ../$MODEL ../$dMODEL ../$DATA wFile_Model wFile_Data -v full]" | tee -a std_out.txt
#
#-A  J rFile_Model rFile_dModel rFile_Data [wFile_Model wFile_Data]
mpirun -n $NCORES ../$EXEC -A J ../$MODEL ../$dMODEL ../$DATA wFile_Model wFile_Data -v full | tee -a std_out.txt
#
# CATCH RESULT
RESULT=$?
#
# TEST RESULT
if [ "$RESULT" -ne "0" ]; then
	#
	#
	echo "TEST ADJ J FAIL: $RESULT" | tee -a std_out.txt
	#
	#
	cd ..
	#
	#
	exit $RESULT
fi
#
#
echo "#### FINISH ADJ J MPI TEST ####" | tee -a std_out.txt
#
#
cd ..
#
#
mv ${EXEC_NAME}/ outputs/temp/test_adj_j
#
#
exit 0
#
# END OF SCRIPT

