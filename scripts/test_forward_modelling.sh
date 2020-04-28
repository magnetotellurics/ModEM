#!/bin/bash
#
# ARGUMENTS: 1 - Mod3DMT PATH, 2 - USER_CTRL FILE, 3 - NUMER OF CORES
#
EXEC_PATH	=$1
USER_CTRL	=$2
ncores		=$3
#
EXEC_NAME="${EXEC_PATH##*/}"
#
# CREATE TEST OUTPUT FOLDER
mkdir -p test_forward_modelling_$EXEC_NAME
#
# ENTER TEST OUTPUT FOLDER
cd test_forward_modelling_$EXEC_NAME/
#
#
echo "#### START FORWARD MODELLING $EXEC_NAME MPI TEST WITH $ncores CORES AT $now ####" | tee std_out.txt
#
#
echo "#### COMMAND LINE: [mpirun -n $ncores ../$EXEC_PATH -W ../$USER_CTRL -v full]" | tee -a std_out.txt
#
#
mpirun -n $ncores ../$EXEC_PATH -W ../$USER_CTRL -v full | tee -a std_out.txt
#
# CATCH RESULT
result=$?
#
# TEST RESULT
if [ "$result" -ne "0" ]; then
	#
	#
	echo "TEST FORWARD MODELLING $EXEC_NAME FAIL: $result" | tee -a std_out.txt
	#
	#
	cd ..
	#
	#
	exit $result
fi
#
#
echo "#### FINISH FORWARD MODELLING $EXEC_PATH MPI TEST ####" | tee -a std_out.txt
#
#
cd ..
#
#
mv test_forward_modelling_$EXEC_NAME/ outputs/temp/
#
#
exit 0
#
# END OF SCRIPT

