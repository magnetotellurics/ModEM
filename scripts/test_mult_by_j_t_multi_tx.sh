#!/bin/bash
#
# ARGUMENTS: 1 - Mod3DMT EXECUTABLE, 2 - MODEL FILE, 3 - DATA FILE, 4 - NUMBER OF CORES

EXEC=$1
MODEL=$2
DATA=$3
ncores=$4
#
# STRING NOW
now=$(date "+%Y/%m/%d - %H:%M:%S")
#
# REMOVE TEST OUTPUT FOLDER FROM MAIN OUTPUT FOLDER
rm -rf outputs/test_mult_by_j_t_multi_tx/
#
# CREATE TEST OUTPUT FOLDER
mkdir -p test_mult_by_j_t_multi_tx
#
# ENTER TEST OUTPUT FOLDER
cd test_mult_by_j_t_multi_tx/
#
#
echo "#### START MULT_BY_J_T_multi_Tx MPI TEST WITH $ncores CORES AT $now ####" | tee std_out.txt
#
#
echo "#### COMMAND LINE: [mpirun -n $ncores ../$EXEC -x ../$MODEL ../$DATA wFile_dModel -v full]" | tee -a std_out.txt
#
#
mpirun -n $ncores ../$EXEC -x ../$MODEL ../$DATA wFile_dModel -v full | tee -a std_out.txt
#
# CATCH RESULT
result=$?
#
# TEST RESULT
if [ "$result" -ne "0" ]; then
	#
	#
	echo "TEST MULT_BY_J_T_multi_Tx FAIL: $result" | tee -a std_out.txt
	#
	#
	cd ..
	#
	#
	exit $result
fi
#
#
echo "#### FINISH MULT_BY_J_T_multi_Tx MPI TEST ####" | tee -a std_out.txt
#
#
cd ..
#
#
mv test_mult_by_j_t_multi_tx/ outputs/
#
#
exit 0
#
# END OF SCRIPT

