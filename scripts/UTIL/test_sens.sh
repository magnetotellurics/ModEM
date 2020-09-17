#!/bin/bash
#
# ARGUMENTS: 1 - Mod3DMT EXECUTABLE, 2 - MODEL FILE, 3 - dMODEL FILE, 4 - DATA FILE, 5 - NUMBER_OF CORES
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
rm -rf outputs/test_sens/
#
# CREATE TEST OUTPUT FOLDER
mkdir -p test_sens
#
# ENTER TEST OUTPUT FOLDER
cd test_sens/
#
#
echo "#### START SENS MPI TEST WITH $ncores CORES AT $now ####" | tee std_out.txt
#
#
echo "#### COMMAND LINE: [mpirun -n $ncores ../$EXEC -S ../$MODEL ../$dMODEL ../$DATA wFile_Data wFile_Sens -v full]" | tee -a std_out.txt
#
#
mpirun -n $ncores ../$EXEC -S ../$MODEL ../$dMODEL ../$DATA wFile_Data wFile_Sens -v full | tee -a std_out.txt
#
# CATCH RESULT
result=$?
#
# TEST RESULT
if [ "$result" -ne "0" ]; then
	#
	#
	echo "TEST SENS FAIL: $result" | tee -a std_out.txt
	#
	#
	cd ..
	#
	#
	exit $result
fi
#
#
echo "#### FINISH SENS MPI TEST ####" | tee -a std_out.txt
#
#
cd ..
#
#
mv test_sens/ outputs/
#
#
exit 0
#
# END OF SCRIPT

