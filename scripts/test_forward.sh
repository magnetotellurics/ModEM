#!/bin/bash
#
# ARGUMENTS: 1 - Mod3DMT EXECUTABLE, 2 - MODEL FILE, 3 - DATA FILE, 4 - NUMER OF CORES
EXEC=$1
MODEL=$2
DATA=$3
ncores=$4
#
# CREATE TEST OUTPUT FOLDER
mkdir -p test_forward_$EXEC
#
# ENTER TEST OUTPUT FOLDER
cd test_forward_$EXEC/
#
#
echo "#### START FORWARD $EXEC MPI TEST WITH $ncores CORES AT $now ####" | tee std_out.txt
#
#
echo "#### COMMAND LINE: [mpirun -n $ncores ../$EXEC -F ../$MODEL ../$DATA wFile_Data.dat wFile_EMsoln -v full]" | tee -a std_out.txt
#
#
mpirun -n $ncores ../$EXEC -F ../$MODEL ../$DATA wFile_Data.dat wFile_EMsoln -v full | tee -a std_out.txt
#
# CATCH RESULT
result=$?
#
# TEST RESULT
if [ "$result" -ne "0" ]; then
	#
	#
	echo "TEST FORWARD $EXEC FAIL: $result" | tee -a std_out.txt
	#
	#
	cd ..
	#
	#
	exit $result
fi
#
#
echo "#### FINISH FORWARD $EXEC MPI TEST ####" | tee -a std_out.txt
#
#
cd ..
#
mv test_forward_$EXEC/ outputs/temp/
#
#
exit 0
#
# END OF SCRIPT

