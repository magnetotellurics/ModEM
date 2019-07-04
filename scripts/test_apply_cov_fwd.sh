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
rm -rf outputs/test_apply_cov_fwd/
#
# CREATE TEST OUTPUT FOLDER
mkdir -p test_apply_cov_fwd
#
# ENTER TEST OUTPUT FOLDER
cd test_apply_cov_fwd/
#
#
echo "#### START TEST_APPLY COV FWD MPI TEST WITH $ncores CORES AT $now ####" | tee std_out.txt
#
#
echo "#### COMMAND LINE: [mpirun -n $ncores ../$EXEC -C FWD ../$MODEL wFile_Model.ws -v full]" | tee std_out.txt
#
#
mpirun -n $ncores ../$EXEC -C FWD ../$MODEL wFile_Model.ws -v full | tee std_out.txt
#
# CATCH RESULT
result=$?
#
# TEST RESULT
if [ "$result" -ne "0" ]; then
	#
	#
	echo "TEST_APPLY COV FWD FAIL: $result" | tee std_out.txt
	#
	#
	cd ..
	#
	#
	exit $result
fi
#
#
echo "#### FINISH TEST_APPLY COV FWD MPI TEST ####" | tee std_out.txt
#
#
cd ..
#
#
mv test_apply_cov_fwd/ outputs/
#
#
exit 0
#
# END OF SCRIPT

