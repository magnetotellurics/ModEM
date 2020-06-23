#!/bin/bash
#
# ARGUMENTS: 1 - Mod3DMT EXECUTABLE, 2 - DATA FILE, 3 - NUMBER OF CORES
EXEC=$1
DATA=$2
ncores=$3
#
# STRING NOW
now=$(date "+%Y/%m/%d - %H:%M:%S")
#
# REMOVE TEST OUTPUT FOLDER FROM MAIN OUTPUT FOLDER
rm -rf outputs/test_adj_d/
#
# CREATE TEST OUTPUT FOLDER
mkdir -p test_adj_d
#
# ENTER TEST OUTPUT FOLDER
cd test_adj_d/
#
#
echo "#### START ADJ D MPI TEST WITH $ncores CORES AT $now ####" | tee std_out.txt
#
#
echo "#### COMMAND LINE: [mpirun -n $ncores ../$EXEC -A d ../$DATA wFile_Data -v full]" | tee -a std_out.txt
#
#-A  d rFile_Data wFile_Data [delta]
mpirun -n $ncores ../$EXEC -A d ../$DATA wFile_Data -v full | tee -a std_out.txt
#
# CATCH RESULT
result=$?
#
# TEST RESULT
if [ "$result" -ne "0" ]; then
	#
	#
	echo "TEST ADJ D FAIL: $result" | tee -a std_out.txt
	#
	#
	cd ..
	#
	#
	exit $result
fi
#
#
echo "#### FINISH ADJ D MPI TEST ####" | tee -a std_out.txt
#
#
cd ..
#
#
mv test_adj_d/ outputs/
#
#
exit 0
#
# END OF SCRIPT

