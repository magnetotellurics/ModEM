#!/bin/bash
#
# ARGUMENTS: 1 - Mod3DMT EXECUTABLE, 2 - MODEL FILE, 3 - EMrhs FILE, 4 - DATA FILE, 5 - NUMBER OF CORES
EXEC=$1
MODEL=$2
EMrhs=$3
DATA=$4
ncores=$5
#
# STRING NOW
now=$(date "+%Y/%m/%d - %H:%M:%S")
#
# REMOVE TEST OUTPUT FOLDER FROM MAIN OUTPUT FOLDER
rm -rf outputs/test_adj_s/
#
# CREATE TEST OUTPUT FOLDER
mkdir -p test_adj_s
#
# ENTER TEST OUTPUT FOLDER
cd test_adj_s/
#
#
echo "#### START ADJ S MPI TEST WITH $ncores CORES AT $now ####" | tee std_out.txt
#
#
echo "#### COMMAND LINE: [mpirun -n $ncores ../$EXEC -A S ../$MODEL ../$EMrhs ../$DATA wFile_EMsoln -v full]" | tee -a std_out.txt
#
#-A  S rFile_Model rFile_EMrhs rFile_Data [wFile_EMsoln]
mpirun -n $ncores ../$EXEC -A S ../$MODEL ../$EMrhs ../$DATA wFile_EMsoln -v full | tee -a std_out.txt
#
# CATCH RESULT
result=$?
#
# TEST RESULT
if [ "$result" -ne "0" ]; then
	#
	#
	echo "TEST ADJ S FAIL: $result" | tee -a std_out.txt
	#
	#
	cd ..
	#
	#
	exit $result
fi
#
#
echo "#### FINISH ADJ S MPI TEST ####" | tee -a std_out.txt
#
#
cd ..
#
#
mv test_adj_s/ outputs/
#
#
exit 0
#
# END OF SCRIPT

