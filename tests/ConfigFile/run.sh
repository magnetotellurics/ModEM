#!/bin/sh
EXEC=$1
echo $1
pwd
ls
ncores=64
echo "COMAND: [mpirun -n $ncores ../$EXEC -W rFile_Config.txt]"
mpirun -n $ncores ../$EXEC -W rFile_Config.txt
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
echo "#### FINISH MPI TEST ####" | tee -a std_out.txt
#.
#
exit 0
#
# END OF SCRIPT

