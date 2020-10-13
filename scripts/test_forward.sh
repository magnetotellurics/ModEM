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
echo "#### START FORWARD MPI TEST WITH $NCORES CORES AT $NOW ####" | tee std_out.txt
#
#
echo "#### COMMAND LINE: [mpirun -n $NCORES ../$EXEC -F ../$MODEL ../$DATA wFile_Data.dat wFile_EMsoln -v full]" | tee -a std_out.txt
#
#
mpirun -n $NCORES ../$EXEC -F ../$MODEL ../$DATA wFile_Data.dat wFile_EMsoln -v full | tee -a std_out.txt
#
# CATCH RESULT
RESULT=$?
#
# TEST RESULT
if [ "$RESULT" -ne "0" ]; then
	#
	#
	echo "TEST FORWARD FAIL: $RESULT" | tee -a std_out.txt
	#
	#
	cd ..
	#
	#
	exit $RESULT
fi
#
#
echo "#### FINISH FORWARD MPI TEST ####" | tee -a std_out.txt
#
#
cd ..
#
#
mv ${EXEC_NAME}/ outputs/temp/test_forward
#
# BUILD bin/SolverDiagnostic3D
cd tools/SolverDiagnostic3D/
bash build_linux.sh
#
cd ../../outputs/temp/test_forward/${EXEC_NAME}
../../../../tools/SolverDiagnostic3D/bin/SolverDiagnostic3D *SolverStatFile_* ../../../../tools/MathBox/mathbox-bundle.js
cd ../../../..
#
#
exit 0
#
# END OF SCRIPT

