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
		mkdir -p outputs/temp/test_symmetry/${DIR_NAME}
		#
		# ENTER TEST OUTPUT FOLDER
		cd outputs/temp/test_symmetry/${DIR_NAME}
		#
		# CREATE TEST OUTPUT FOLDER
		mkdir -p ${EXEC_NAME}
		#
		# ENTER TEST OUTPUT FOLDER
		cd ${EXEC_NAME}
		#
		T_START=$(date +%s%3N)
		#
		echo "	> [${DIR_NAME}]: START TEST SYMMETRY ${EXEC_NAME} WITH $NCORES CORES" | tee -a ../../../summary.txt
		#
		echo "#### [${DIR_NAME}]:START SYMMETRY ${EXEC_NAME} TEST WITH $NCORES CORES AT $NOW ####" | tee std_out.txt
		#
		echo "#### COMMAND LINE: [mpirun -n $NCORES ../../../../../$EXEC -A J ../../../../../$DIR/rFile_Model ../../../test_mult_by_j_t/${DIR_NAME}/${EXEC_NAME}/wFile_dModel ../../../../../$DIR/rFile_Data -v full]" | tee -a std_out.txt
		#
		mpirun -n $NCORES ../../../../../$EXEC -A J ../../../../../$DIR/rFile_Model ../../../test_mult_by_j_t/${DIR_NAME}/${EXEC_NAME}/wFile_dModel ../../../../../$DIR/rFile_Data -v full | tee -a std_out.txt
		#
		# CATCH RESULT
		RESULT=$?
		#
		# TEST RESULT
		if [ "$RESULT" -ne "0" ]; then
			echo "	> [${DIR_NAME}]: TEST ADJ J ${EXEC_NAME} FAIL: $RESULT" | tee -a ../../../summary.txt
		else
			echo "	> [${DIR_NAME}]: TEST ADJ J ${EXEC_NAME} PASS: $RESULT" | tee -a ../../../summary.txt
		fi
		#
		T_END=$(date +%s%3N)
		echo "	> [${DIR_NAME}]: Time Spent: $(( ( $T_END - $T_START ) / 1000 )) seconds" | tee -a ../../../summary.txt
		echo "	#" | tee -a ../../../summary.txt
		#
		# CREATE FILE WITH 2 dot product LINES
		grep -hnr "dot product" std_out.txt > dot_product.txt
		#
		# GET VALUES
		DOT1=$(sed '1q;d' dot_product.txt | awk '{ print $5 }')
		DOT2=$(sed '2q;d' dot_product.txt | awk '{ print $5 }')
		#
		D1=$(echo $DOT1 | cut -c1-8)
		D2=$(echo $DOT2 | cut -c1-8)
		#
		# SYMMETRY TEST
		if [[ -n "$D1" && -n "$D2" && "$D1" == "$D2" ]]; then
			echo "	> [${DIR_NAME}]: TEST SYMMETRY ${EXEC_NAME} PASS: $D1 == $D2" | tee -a ../../../summary.txt
		else
			echo "	> [${DIR_NAME}]: TEST SYMMETRY ${EXEC_NAME} FAIL: $D1 != $D2" | tee -a ../../../summary.txt
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
exit 0
#
# END OF SCRIPT

