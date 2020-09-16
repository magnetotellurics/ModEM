#!/bin/bash
#
# ARGUMENTS: 1 - Mod3DMT EXECUTABLE, 2 - MODEL FILE, 3 - DATA FILE, 4- EMrhs FILE, 5 - NUMBER OF CORES
EXEC=$1
MODEL=$2
DATA=$3
EMrhs=$4
ncores=$5
#
# STRING NOW
now=$(date "+%Y/%m/%d - %H:%M:%S")

#
# REMOVE TEST OUTPUT FOLDER FROM MAIN OUTPUT FOLDER
rm -rf outputs/test_adj_b/
#
# CREATE TEST OUTPUT FOLDER
mkdir -p test_adj_b
#
# ENTER TEST OUTPUT FOLDER
cd test_adj_b/
#
#
echo "#### START ADJ b MPI TEST WITH $ncores CORES AT $now ####" | tee std_out.txt
#
#
echo "#### COMMAND LINE: [mpirun -n $ncores ../$EXEC -A b ../$MODEL ../$DATA ../$EMrhs  wFile_EMrhs -v full]" | tee -a std_out.txt
#
#-A  b rFile_Model rFile_Data rFile_EMrhs wFile_EMrhs [delta]
mpirun -n $ncores ../$EXEC -A b ../$MODEL ../$DATA ../$EMrhs  wFile_EMrhs -v full | tee -a std_out.txt
#
# CATCH RESULT
result=$?
#
# TEST RESULT
if [ "$result" -ne "0" ]; then
	#
	#
	echo "TEST ADJ b FAIL: $result" | tee -a std_out.txt
	#
	#
	cd ..
	#
	#
	exit $result
fi
#
#
echo "#### FINISH ADJ b MPI TEST ####" | tee -a std_out.txt
#
#
cd ..
#
#
mv test_adj_b/ outputs/temp/
#
#
exit 0
#
# END OF SCRIPT

