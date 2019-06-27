#!/bin/bash
#
# ARGUMENTS: 1 - Mod3DMT EXECUTABLE, 2 - MODEL FILE
EXEC=$1
MODEL=$2
#
# STRING NOW
now=$(date "+%Y/%m/%d - %H:%M:%S")
#
# GET ENVIROMENT NUMBER OF CORES 
ncores=$(nproc)
#
# REMOVE TEST OUTPUT FOLDER FROM MAIN OUTPUT FOLDER
rm -rf outputs/test_cov/
#
# CREATE TEST OUTPUT FOLDER
mkdir -p test_cov
#
# ENTER TEST OUTPUT FOLDER
cd test_cov/
#
#
echo "#### START TEST_COV MPI TEST WITH $ncores CORES AT $now ####\n" >> std_out.txt
#
#
echo "#### COMMAND LINE: [mpirun -n $ncores ../$EXEC -C ../$MODEL wFile_Model.ws -v full]" >> std_out.txt
#
#
mpirun -n $ncores ../$EXEC -C ../$MODEL wFile_Model.ws -v full &>> std_out.txt
#
# CATCH RESULT
result=$?
#
# TEST RESULT
if [ "$result" -ne "0" ]; then
	#
	#
	echo "TEST_COV FAIL: $result"
	#
	#
	cd ..
	#
	#
	exit $result
fi
#
#
echo "### FINISH TEST_COV MPI TEST ###\n" >> std_out.txt
#
#
cd ..
#
#
mv test_cov/ outputs/
#
#
# END OF SCRIPT

