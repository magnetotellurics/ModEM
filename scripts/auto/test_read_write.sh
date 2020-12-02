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
		# CREATE AND ENTER TEST OUTPUT FOLDER
		mkdir -p outputs/temp/test_read_write/${DIR_NAME}/${EXEC_NAME} && cd "$_"
		#
		echo $(pwd)
		#
		T_START=$(date +%s%3N)
		#
		echo "#### [${DIR_NAME}]: START TEST READ_WRITE ${EXEC_NAME} WITH $NCORES CORES AT $NOW ####" | tee -a std_out.txt
		#
		echo "	> [${DIR_NAME}]: START TEST READ_WRITE ${EXEC_NAME} WITH $NCORES CORES" | tee -a ../../../summary.txt
		#
		echo "#### COMMAND LINE: [mpirun -n $NCORES ../../../../../$EXEC -R ../../../../../$DIR/rFile_Model ../../../../../$DIR/rFile_Data wFile_Model wFile_Data -v full]" | tee -a std_out.txt
		#
		mpirun -n $NCORES ../../../../../$EXEC -R ../../../../../$DIR/rFile_Model ../../../../../$DIR/rFile_Data wFile_Model wFile_Data -v full | tee -a std_out.txt
		#
		# CATCH RESULT
		RESULT=$?
		#
		# TEST RESULT
		if [ "$RESULT" -ne "0" ]; then
			#
			echo "	> [${DIR_NAME}]: TEST READ_WRITE ${EXEC_NAME} FAIL: $RESULT" | tee -a ../../../summary.txt
			T_END=$(date +%s%3N)
			echo "	> [${DIR_NAME}]: Time Spent: $(( ( $T_END - $T_START ) / 1000 )) seconds" | tee -a ../../../summary.txt
			echo "	#" | tee -a ../../../summary.txt
			#
			cd ../../../../..
			#
			exit $RESULT
		fi
		#
		echo "#### [${DIR_NAME}]: FINISH READ_WRITE ${EXEC_NAME} TEST ####" | tee -a std_out.txt
		#
		echo "	> [${DIR_NAME}]: TEST READ_WRITE ${EXEC_NAME} PASS: $RESULT" | tee -a ../../../summary.txt
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

