#!/bin/bash
#
# ARGUMENTS: 1 - Mod3DMT EXECUTABLE, 2 - MODEL FILE, 3 - DATA FILE, 4- EMsoln FILE
EXEC=$1
MODEL=$2
DATA=$3
EMsoln=$4
#
# STRING NOW
now=$(date "+%Y/%m/%d - %H:%M:%S")
#
# GET ENVIROMENT NUMBER OF CORES 
ncores=$(nproc)
#
# REMOVE TEST OUTPUT FOLDER FROM MAIN OUTPUT FOLDER
rm -rf outputs/test_adj_e/
#
# CREATE TEST OUTPUT FOLDER
mkdir -p test_adj_e
#
# ENTER TEST OUTPUT FOLDER
cd test_adj_e/
#
#
echo "#### START ADJ e MPI TEST WITH $ncores CORES AT $now ####" | tee std_out.txt
#
#
echo "#### COMMAND LINE: [mpirun -n $ncores ../$EXEC -A e ../$MODEL ../$DATA ../$EMsoln wFile_EMsoln -v full]" | tee std_out.txt
#
#-A  e rFile_Model rFile_Data rFile_EMsoln wFile_EMsoln [delta]
mpirun -n $ncores ../$EXEC -A e ../$MODEL ../$DATA ../$EMsoln  wFile_EMsoln -v full | tee std_out.txt
#
# CATCH RESULT
result=$?
#
# TEST RESULT
if [ "$result" -ne "0" ]; then
	#
	#
	echo "TEST ADJ e FAIL: $result" | tee std_out.txt
	#
	#
	cd ..
	#
	#
	exit $result
fi
#
#
echo "#### FINISH ADJ e MPI TEST ####" | tee std_out.txt
#
#
cd ..
#
#
mv test_adj_e/ outputs/
#
#
exit 0
#
# END OF SCRIPT

