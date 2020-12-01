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
emptyspace=""
#
# GET ENVIROMENT NUMBER OF CORES 
NCORES=$(nproc)
#
DIRS=inputs/*
for DIR in $DIRS; do
	#
    if [ -d "$DIR" ]; then
		#
		DIR_NAME=$DIR
		DIR_NAME=${DIR_NAME/*\//$emptyspace}
		#
        echo $DIR_NAME
		#
		# CREATE TEST OUTPUT FOLDER
		mkdir -p ${DIR_NAME}
		#
		# ENTER TEST OUTPUT FOLDER
		cd ${DIR_NAME}/
		#
		# CREATE TEST OUTPUT FOLDER
		mkdir -p ${EXEC_NAME}
		#
		# ENTER TEST OUTPUT FOLDER
		cd ${EXEC_NAME}/
		#
		echo "#### START MODEM MPI TEST WITH $NCORES CORES AT $NOW ####" | tee -a std_out.txt
		#
		echo "#### COMMAND LINE: [mpirun -n $NCORES ../../$EXEC]" | tee -a std_out.txt
		#
		mpirun -n $NCORES ../../$EXEC | tee -a std_out.txt
		#
		# CATCH RESULT
		RESULT=$?
		#
		# TEST RESULT
		if [ "$RESULT" -ne "0" ]; then
			#
			echo "TEST PLAIN FAIL: $RESULT" | tee -a std_out.txt
			#
			cd ../..
			#
			exit $RESULT
		fi
		#
		echo "#### FINISH PAIN MPI TEST ####" | tee -a std_out.txt
		#
		cd ../..
		#
		mv ${DIR_NAME}/ outputs/temp/test_plain/
    fi
done
#
#
exit 0
#
# END OF SCRIPT

