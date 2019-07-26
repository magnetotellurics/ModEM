#!/bin/sh
EXEC=$1
echo $1
pwd
ncores=64
mpirun -n $ncores ../$EXEC -W rFile_config.txt
#
# CATCH RESULT
result=$?
#
# TEST RESULT
if [ "$result" -ne "0" ]; then
	#
	#
	echo "TEST FAIL: $result" | tee -a std_out.txt
	#
	exit $result
fi
#
#
echo "#### FINISH GRAD MPI TEST ####" | tee -a std_out.txt
#.
#
exit 0
#
# END OF SCRIPT

