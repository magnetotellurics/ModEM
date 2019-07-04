#!/bin/bash
#
# ARGUMENTS: 1 - Mod3DMT EXECUTABLE, 2 - dMODEL FILE, 3 - DATA FILE, 4 - DATA FILE, 5 - NUMBER OF CORES
EXEC=$1
MODEL=$2
dMODEL=$3
DATA=$4
ncores=$5
#
# STRING NOW
now=$(date "+%Y/%m/%d - %H:%M:%S")
#
# REMOVE TEST OUTPUT FOLDER FROM MAIN OUTPUT FOLDER
rm -rf outputs/test_mult_by_j/
#
# CREATE TEST OUTPUT FOLDER
mkdir -p test_mult_by_j
#
# ENTER TEST OUTPUT FOLDER
cd test_mult_by_j/
#
#
echo "#### START MULT_BY_J MPI TEST WITH $ncores CORES AT $now ####" | tee std_out.txt
#
#
echo "#### COMMAND LINE: [mpirun -n $ncores ../$EXEC -M ../$MODEL ../$dMODEL ../$DATA wFile_Data -v full]" | tee std_out.txt
#
#
mpirun -n $ncores ../$EXEC -M ../$MODEL ../$dMODEL ../$DATA wFile_Data -v full | tee std_out.txt
#
# CATCH RESULT
result=$?
#
# TEST RESULT
if [ "$result" -ne "0" ]; then
	#
	#
	echo "TEST MULT_BY_J FAIL: $result" | tee std_out.txt
	#
	#
	cd ..
	#
	#
	exit $result
fi
#
#
echo "#### FINISH MULT_BY_J MPI TEST ####" | tee std_out.txt
#
#
cd ..
#
#
mv test_mult_by_j/ outputs/
#
#
exit 0
#
# END OF SCRIPT

