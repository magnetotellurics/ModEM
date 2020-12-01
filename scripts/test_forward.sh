#!/bin/bash
#
# ARGUMENTS: 1 - Mod3DMT EXECUTABLE
EXEC=$1
EXEC_NAME=$EXEC
emptyspace=""
EXEC_NAME=${EXEC_NAME/.exe/$emptyspace}
EXEC_NAME=${EXEC_NAME/.sh/$emptyspace}
EXEC_NAME=${EXEC_NAME/.txt/$emptyspace}
EXEC_NAME=${EXEC_NAME/*\//$emptyspace}
#
#
if [[ $2 -eq 0 ]]; then
	NCORES=$(nproc)
else
	NCORES=$2
fi
#
emptyspace=""
#
# STRING NOW
NOW=$(date "+%Y/%m/%d - %H:%M:%S")
#
DIRS=inputs/*
for DIR in $DIRS; do
	#
    if [ -d "$DIR" ]; then
		#
		DIR_NAME=$DIR
		DIR_NAME=${DIR_NAME/*\//$emptyspace}
		#
		# CREATE TEST OUTPUT FOLDER
		mkdir -p outputs/temp/test_forward/${DIR_NAME}
		#
		# ENTER TEST OUTPUT FOLDER
		cd outputs/temp/test_forward/${DIR_NAME}
		#
		# CREATE TEST OUTPUT FOLDER
		mkdir -p ${EXEC_NAME}
		#
		# ENTER TEST OUTPUT FOLDER
		cd ${EXEC_NAME}
		#
		T_START=$(date +%s%3N)
		#
		echo "	> [${DIR_NAME}]: START TEST FORWARD ${EXEC_NAME} WITH $NCORES CORES" | tee -a ../../../summary.txt
		#
		echo "#### [${DIR_NAME}]:START FORWARD ${EXEC_NAME} TEST WITH $NCORES CORES AT $NOW ####" | tee -a std_out.txt
		#
		echo "#### COMMAND LINE: [mpirun -n $NCORES ../../$EXEC -F ../../outputs/temp/test_read_write/${DIR_NAME}/${EXEC_NAME}/wFile_Model ../../outputs/temp/test_read_write/${DIR_NAME}/${EXEC_NAME}/wFile_Data wFile_Data wFile_EMsoln -v full]" | tee -a std_out.txt
		#
		mpirun -n $NCORES ../../../../../$EXEC -F ../../../../../$DIR/rFile_Model ../../../../../$DIR/rFile_Data wFile_Data wFile_EMsoln -v full | tee -a std_out.txt
		#
		# CATCH RESULT
		RESULT=$?
		#
		# TEST RESULT
		if [ "$RESULT" -ne "0" ]; then
			#
			echo "	> [${DIR_NAME}]: TEST FORWARD ${EXEC_NAME} FAIL: $RESULT" | tee -a ../../../summary.txt
			T_END=$(date +%s%3N)
			echo "	> [${DIR_NAME}]: Time Spent: $(( ( $T_END - $T_START ) / 1000 )) seconds" | tee -a ../../../summary.txt
			echo "	#" | tee -a ../../../summary.txt
			#
			cd ../../../../..
			#
			exit $RESULT
		fi
		#
		echo "	> [${DIR_NAME}]: TEST FORWARD ${EXEC_NAME} PASS: $RESULT" | tee -a ../../../summary.txt
		T_END=$(date +%s%3N)
		echo "	> [${DIR_NAME}]: Time Spent: $(( ( $T_END - $T_START ) / 1000 )) seconds" | tee -a ../../../summary.txt
		#
		echo "	> [${DIR_NAME}]: CREATE SolverDiagnostic3D" | tee -a ../../../summary.txt
		#
		T_START=$(date +%s%3N)
		#
		cd ../../../../..
		#
		# BUILD bin/SolverDiagnostic3D
		cd tools/SolverDiagnostic3D/
		bash build_linux.sh
		#
		cd ../../outputs/temp/test_forward/${DIR_NAME}/${EXEC_NAME}
		../../../../../tools/SolverDiagnostic3D/bin/SolverDiagnostic3D *_SolverStatFile_* ../../../../../tools/MathBox/mathbox-bundle.js
		rm -rf ../../../../../tools/SolverDiagnostic3D/bin
		cd ../../../../..
		#
		T_END=$(date +%s%3N)
		echo "	> [${DIR_NAME}]: Time Spent: $(( ( $T_END - $T_START ) / 1000 )) seconds" | tee -a outputs/temp/summary.txt
		echo "	#"  | tee -a outputs/temp/summary.txt
	fi
done
#
exit 0
#
# END OF SCRIPT

