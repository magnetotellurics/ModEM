#!/bin/bash
#
# ARGUMENTS: 1 - Mod3DMT EXECUTABLE, 2 - MODEL FILE, 3 - dMODEL FILE, 4 - DATA FILE, 5 - NUMBER OF CORES
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
rm -rf outputs/test_adj_q/
#
# CREATE TEST OUTPUT FOLDER
mkdir -p test_adj_q
#
# ENTER TEST OUTPUT FOLDER
cd test_adj_q/
#
#
echo "#### START ADJ Q MPI TEST WITH $ncores CORES AT $now ####" | tee std_out.txt
#
#
echo "#### COMMAND LINE: [mpirun -n $ncores ../$EXEC -A Q ../$MODEL ../$dMODEL ../$DATA wFile_Model wFile_Data -v full]" | tee -a std_out.txt
#
#-A  Q rFile_Model rFile_dModel rFile_Data [wFile_Model wFile_Data]
mpirun -n $ncores ../$EXEC -A Q ../$MODEL ../$dMODEL ../$DATA wFile_Model wFile_Data -v full | tee -a std_out.txt
#
# CATCH RESULT
result=$?
#
# TEST RESULT
if [ "$result" -ne "0" ]; then
	#
	#
	echo "TEST ADJ Q FAIL: $result" | tee -a std_out.txt
	#
	#
	cd ..
	#
	#
	exit $result
fi
#
#
echo "#### FINISH ADJ Q MPI TEST ####" | tee -a std_out.txt
#
#
cd ..
#
#
mv test_adj_q/ outputs/
#
#
exit 0
#
# END OF SCRIPT

