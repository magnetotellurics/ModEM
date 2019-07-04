#!/bin/bash
#
# ARGUMENTS: 1 - Mod3DMT EXECUTABLE, 2 - MODEL FILE, 3 - NUMBER OF CORES
EXEC=$1
MODEL=$2
ncores=$3
#
# STRING NOW
now=$(date "+%Y/%m/%d - %H:%M:%S")
#
# REMOVE TEST OUTPUT FOLDER FROM MAIN OUTPUT FOLDER
rm -rf outputs/test_adj_m/
#
# CREATE TEST OUTPUT FOLDER
mkdir -p test_adj_m
#
# ENTER TEST OUTPUT FOLDER
cd test_adj_m/
#
#
echo "#### START ADJ M MPI TEST WITH $ncores CORES AT $now ####" | tee std_out.txt
#
#
echo "#### COMMAND LINE: [mpirun -n $ncores ../$EXEC -A m ../$MODEL wFile_Model -v full]" | tee std_out.txt
#
#-A  m rFile_Model wFile_Model [delta]
mpirun -n $ncores ../$EXEC -A m ../$MODEL wFile_Model -v full | tee std_out.txt
#
# CATCH RESULT
result=$?
#
# TEST RESULT
if [ "$result" -ne "0" ]; then
	#
	#
	echo "TEST ADJ M FAIL: $result" | tee std_out.txt
	#
	#
	cd ..
	#
	#
	exit $result
fi
#
#
echo "#### FINISH ADJ M MPI TEST ####" | tee std_out.txt
#
#
cd ..
#
#
mv test_adj_m/ outputs/
#
#
exit 0
#
# END OF SCRIPT

