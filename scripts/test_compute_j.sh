#!/bin/bash
#
# ARGUMENTS: 1 - Mod3DMT EXECUTABLE, 2 - MODEL FILE, 3 - DATA FILE, 4 - NUMBER OF CORES
EXEC=$1
MODEL=$2
DATA=$3
ncores=$4
#
# STRING NOW
now=$(date "+%Y/%m/%d - %H:%M:%S")
#
# REMOVE TEST OUTPUT FOLDER FROM MAIN OUTPUT FOLDER
rm -rf outputs/test_compute_j/
#
# CREATE TEST OUTPUT FOLDER
mkdir -p test_compute_j
#
# ENTER TEST OUTPUT FOLDER
cd test_compute_j/
#
#
echo "#### START COMPUT_J MPI TEST WITH $ncores CORES AT $now ####" | tee std_out.txt
#
#
echo "#### COMMAND LINE: [mpirun -n $ncores ../$EXEC -J ../$MODEL ../$DATA wFile_Sens -v full]" | tee std_out.txt
#
#
mpirun -n $ncores ../$EXEC -J ../$MODEL ../$DATA wFile_Sens -v full | tee std_out.txt
#
# CATCH RESULT
result=$?
#
# TEST RESULT
if [ "$result" -ne "0" ]; then
	#
	#
	echo "TEST COMPUT_J FAIL: $result" | tee std_out.txt
	#
	#
	cd ..
	#
	#
	exit $result
fi
#
#
echo "#### FINISH COMPUT_J MPI TEST ####" | tee std_out.txt
#
#
cd ..
#
#
mv test_compute_j/ outputs/
#
#
exit 0
#
# END OF SCRIPT

