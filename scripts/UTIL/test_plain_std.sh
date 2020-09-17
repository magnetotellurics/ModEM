#!/bin/bash
#
# ARGUMENTS: 1 - NUMBER OF CORES
#
# GET ENVIROMENT NUMBER OF CORES 
ncores=$1
#
# STRING NOW
now=$(date "+%Y/%m/%d - %H:%M:%S")
#
#
echo "#### START STD PLAIN MPI TEST WITH $ncores CORES AT $now ####"
#
#
mpirun -n $ncores ./bin/Mod3DMT_STD
#
# CATCH RESULT
result=$?
#
# TEST RESULT
if [ "$result" -ne "0" ]; then
	#
	#
	echo "TEST STD PLAIN FAIL: $result"
	#
	#
	exit $result
fi
#
echo "#### FINISH STD PLAIN TEST ####"
#
#
exit 0
#
# END OF SCRIPT

