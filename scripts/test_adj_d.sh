#!/bin/bash
#
# ARGUMENTS: 1 - Mod3DMT EXECUTABLE, 2 - DATA FILE
EXEC=$1
DATA=$2
#
# STRING NOW
now=$(date "+%Y/%m/%d - %H:%M:%S")
#
# GET ENVIROMENT NUMBER OF CORES 
ncores=$(nproc)
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
echo "#### START ADJ D MPI TEST WITH $ncores CORES AT $now ####" >> std_out.txt
#
#
echo "#### COMMAND LINE: [mpirun -n $ncores ../$EXEC -A d ../$DATA wFile_Data -v full]" >> std_out.txt
#
#-A  d rFile_Data wFile_Data [delta]
mpirun -n $ncores ../$EXEC -A d ../$DATA wFile_Data -v full &>> std_out.txt
#
# CATCH RESULT
result=$?
#
# TEST RESULT
if [ "$result" -ne "0" ]; then
	#
	#
	echo "TEST ADJ D FAIL: $result"
	#
	#
	cd ..
	#
	#
	exit $result
fi
#
#
echo "#### FINISH ADJ D MPI TEST ####" >> std_out.txt
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

