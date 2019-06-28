#!/bin/bash
#
# ARGUMENTS: 1 - Mod3DMT EXECUTABLE, 2 - MODEL FILE, 3 - dMODEL FILE, 4 - EMsoln FILE, 5 - DATA FILE
EXEC=$1
MODEL=$2
dMODEL=$3
EMsoln=$4
DATA=$5
#
# STRING NOW
now=$(date "+%Y/%m/%d - %H:%M:%S")
#
# GET ENVIROMENT NUMBER OF CORES 
ncores=$(nproc)
#
# REMOVE TEST OUTPUT FOLDER FROM MAIN OUTPUT FOLDER
rm -rf outputs/test_adj_p/
#
# CREATE TEST OUTPUT FOLDER
mkdir -p test_adj_p
#
# ENTER TEST OUTPUT FOLDER
cd test_adj_p/
#
#
echo "#### START ADJ P MPI TEST WITH $ncores CORES AT $now ####" >> std_out.txt
#
#
echo "#### COMMAND LINE: [mpirun -n $ncores ../$EXEC -A P ../$MODEL ../$dMODEL ../$EMsoln ../$DATA wFile_Model wFile_EMrhs -v full]" >> std_out.txt
#
#-A  P rFile_Model rFile_dModel rFile_EMsoln rFile_Data [wFile_Model wFile_EMrhs]
mpirun -n $ncores ../$EXEC -A P ../$MODEL ../$dMODEL ../$EMsoln ../$DATA wFile_Model wFile_EMrhs -v full &>> std_out.txt
#
# CATCH RESULT
result=$?
#
# TEST RESULT
if [ "$result" -ne "0" ]; then
	#
	#
	echo "TEST ADJ P FAIL: $result"
	#
	#
	cd ..
	#
	#
	exit $result
fi
#
#
echo "#### FINISH ADJ P MPI TEST ####" >> std_out.txt
#
#
cd ..
#
#
mv test_adj_p/ outputs/
#
#
exit 0
#
# END OF SCRIPT

