#!/bin/bash
#
# ARGUMENTS: 1 - Mod3DMT EXECUTABLE, 2 - MODEL FILE, 3 - dMODEL FILE, 4 - DATA FILE
EXEC=$1
MODEL=$2
dMODEL=$3
DATA=$4
#
# STRING NOW
now=$(date "+%Y/%m/%d - %H:%M:%S")
#
# GET ENVIROMENT NUMBER OF CORES 
ncores=$(nproc)
#
# REMOVE TEST OUTPUT FOLDER FROM MAIN OUTPUT FOLDER
rm -rf outputs/test_adj_j/
#
# CREATE TEST OUTPUT FOLDER
mkdir -p test_adj_j
#
# ENTER TEST OUTPUT FOLDER
cd test_adj_j/
#
#
echo "#### START ADJ J MPI TEST WITH $ncores CORES AT $now ####" >> std_out.txt
#
#
echo "#### COMMAND LINE: [mpirun -n $ncores ../$EXEC -A J ../$MODEL ../$dMODEL ../$DATA wFile_Model wFile_Data -v full]" >> std_out.txt
#
#-A  J rFile_Model rFile_dModel rFile_Data [wFile_Model wFile_Data]
mpirun -n $ncores ../$EXEC -A J ../$MODEL ../$dMODEL ../$DATA wFile_Model wFile_Data -v full &>> std_out.txt
#
# CATCH RESULT
result=$?
#
# TEST RESULT
if [ "$result" -ne "0" ]; then
	#
	#
	echo "TEST ADJ J FAIL: $result"
	#
	#
	cd ..
	#
	#
	exit $result
fi
#
#
echo "#### FINISH ADJ J MPI TEST ####" >> std_out.txt
#
#
cd ..
#
#
mv test_adj_j/ outputs/
#
#
exit 0
#
# END OF SCRIPT

