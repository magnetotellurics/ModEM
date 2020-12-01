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
		mkdir -p outputs/temp/test_mult_by_j_t/${DIR_NAME}
		#
		# ENTER TEST OUTPUT FOLDER
		cd outputs/temp/test_mult_by_j_t/${DIR_NAME}
		#
		# CREATE TEST OUTPUT FOLDER
		mkdir -p ${EXEC_NAME}
		#
		# ENTER TEST OUTPUT FOLDER
		cd ${EXEC_NAME}
		#
		T_START=$(date +%s%3N)
		#
		echo "	> [${DIR_NAME}]: START TEST MULT_BY_J_T ${EXEC_NAME} WITH $NCORES CORES" | tee -a ../../../summary.txt
		#
		echo "#### [${DIR_NAME}]:START MULT_BY_J_T ${EXEC_NAME} TEST WITH $NCORES CORES AT $NOW ####" | tee std_out.txt
		#
		echo "#### COMMAND LINE: [mpirun -n $NCORES ../../../../../$EXEC -T ../../../../../$DIR/rFile_Model ../../../../../$DIR/rFile_Data wFile_dModel -v full]" | tee -a std_out.txt
		#
		mpirun -n $NCORES ../../../../../$EXEC -T ../../../../../$DIR/rFile_Model ../../../../../$DIR/rFile_Data wFile_dModel -v full | tee -a std_out.txt
		#
		# CATCH RESULT
		RESULT=$?
		#
		DMODEL=wFile_dModel
		#
		# TEST RESULT
		if [[ ! -f "$DMODEL" || "$RESULT" -ne "0" ]]; then
			echo "	> [${DIR_NAME}]: TEST MULT_BY_J_T ${EXEC_NAME} FAIL: $RESULT" | tee -a ../../../summary.txt
		else
			echo "	> [${DIR_NAME}]: TEST MULT_BY_J_T ${EXEC_NAME} PASS: $RESULT" | tee -a ../../../summary.txt
		fi
		#
		T_END=$(date +%s%3N)
		echo "	> [${DIR_NAME}]: Time Spent: $(( ( $T_END - $T_START ) / 1000 )) seconds" | tee -a ../../../summary.txt
		echo "	#" | tee -a ../../../summary.txt
		#
		cd ../../../../..
	fi
done
#
#
exit 0
#
# END OF SCRIPT

