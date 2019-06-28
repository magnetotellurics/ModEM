#!/bin/bash
#
# NO ARGUMENTS: 
#
# STRING NOW
now=$(date "+%Y/%m/%d - %H:%M:%S")
#
# GET ENVIROMENT NUMBER OF CORES 
ncores=$(nproc)
#
#
echo "#### START PLAIN MPI TEST WITH $ncores CORES AT $now ####"
#
#
mpirun -n $ncores ./src/Mod3DMT
#
# CATCH RESULT
result=$?
#
# TEST RESULT
if [ "$result" -ne "0" ]; then
	#
	#
	echo "TEST PLAIN FAIL: $result"
	#
	#
	exit $result
fi
#
echo "#### FINISH PLAIN TEST ####"
#
#
exit 0
#
# END OF SCRIPT

