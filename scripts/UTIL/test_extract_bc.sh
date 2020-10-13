#!/bin/bash
#
# ARGUMENTS: 1 - Mod3DMT EXECUTABLE, 2 - MODEL FILE, 3 - DATA FILE, 4 - NUMER OF CORES
EXEC=$1
MODEL=$2
DATA=$3
ncores=$4
#
# STRING NOW
now=$(date "+%Y/%m/%d - %H:%M:%S")
#
# REMOVE TEST OUTPUT FOLDER FROM MAIN OUTPUT FOLDER
rm -rf outputs/test_extract_bc/
#
# CREATE TEST OUTPUT FOLDER
mkdir -p test_extract_bc
#
# ENTER TEST OUTPUT FOLDER
cd test_extract_bc/
#
#
echo "#### START EXTRACT BC MPI TEST WITH $ncores CORES AT $now ####" | tee std_out.txt
#
#
echo "#### COMMAND LINE: [mpirun -n $ncores ../$EXEC -b ../$MODEL ../$DATA wFile_EMrhs -v full]" | tee -a std_out.txt
#
#
mpirun -n $ncores ../$EXEC -b ../$MODEL ../$DATA wFile_EMrhs -v full | tee -a std_out.txt
#
# CATCH RESULT
result=$?
#
# TEST RESULT
if [ "$result" -ne "0" ]; then
	#
	#
	echo "TEST EXTRACT BC FAIL: $result" | tee -a std_out.txt
	#
	#
	cd ..
	#
	#
	exit $result
fi
#
#
echo "#### FINISH EXTRACT BC MPI TEST ####" | tee -a std_out.txt
#
#
cd ..
#
#
mv test_extract_bc/ outputs/
#
#
exit 0
#
# END OF SCRIPT

