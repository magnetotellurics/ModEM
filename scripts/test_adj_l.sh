#!/bin/bash
#
# ARGUMENTS: 1 - Mod3DMT EXECUTABLE, 2 - MODEL FILE, 3 - EMsoln FILE, 3 - DATA FILE
EXEC=$1
MODEL=$2
EMsoln=$3
DATA=$4
#
# STRING NOW
now=$(date "+%Y/%m/%d - %H:%M:%S")
#
# GET ENVIROMENT NUMBER OF CORES 
ncores=$(nproc)
#
# REMOVE TEST OUTPUT FOLDER FROM MAIN OUTPUT FOLDER
rm -rf outputs/test_adj_l/
#
# CREATE TEST OUTPUT FOLDER
mkdir -p test_adj_l
#
# ENTER TEST OUTPUT FOLDER
cd test_adj_l/
#
#
echo "#### START ADJ L MPI TEST WITH $ncores CORES AT $now ####" | tee std_out.txt
#
#
echo "#### COMMAND LINE: [mpirun -n $ncores ../$EXEC -A L ../$MODEL ../$EMsoln ../$DATA wFile_EMrhs wFile_Data -v full]" | tee std_out.txt
#
#-A  L rFile_Model rFile_EMsoln rFile_Data [wFile_EMrhs wFile_Data]
mpirun -n $ncores ../$EXEC -A L ../$MODEL ../$EMsoln ../$DATA wFile_EMrhs wFile_Data -v full | tee std_out.txt
#
# CATCH RESULT
result=$?
#
# TEST RESULT
if [ "$result" -ne "0" ]; then
	#
	#
	echo "TEST ADJ L FAIL: $result" | tee std_out.txt
	#
	#
	cd ..
	#
	#
	exit $result
fi
#
#
echo "#### FINISH ADJ L MPI TEST ####" | tee std_out.txt
#
#
cd ..
#
#
mv test_adj_l/ outputs/
#
#
exit 0
#
# END OF SCRIPT

