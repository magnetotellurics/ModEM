#!/bin/bash
#
# ARGUMENTS: 1 - Mod3DMT EXECUTABLE, 2 - MODEL FILE, 3 - DATA FILE
EXEC=$1
MODEL=$2
DATA=$3
#
# STRING NOW
now=$(date "+%Y/%m/%d - %H:%M:%S")
#
# GET ENVIROMENT NUMBER OF CORES 
ncores=$(nproc)
#
# REMOVE TEST OUTPUT FOLDER FROM MAIN OUTPUT FOLDER
rm -rf outputs/test_inverse/
#
# CREATE TEST OUTPUT FOLDER
mkdir -p test_inverse
#
# ENTER TEST OUTPUT FOLDER
cd test_inverse/
#
#
echo "#### START INVERSE MPI TEST WITH $ncores CORES AT $now ####" | tee std_out.txt
#
#
echo "#### COMMAND LINE: [mpirun -n $ncores ../$EXEC -I NLCG ../$MODEL ../$DATA -v full]" | tee std_out.txt
#
#
mpirun -n $ncores ../$EXEC -I NLCG ../$MODEL ../$DATA -v full | tee std_out.txt
#
# CATCH RESULT
result=$?
#
# TEST RESULT
if [ "$result" -ne "0" ]; then
	#
	#
	echo "TEST FORWARD FAIL: $result" | tee std_out.txt
	#
	#
	cd ..
	#
	#
	exit $result
fi
#
#
echo "#### FINISH FORWARD MPI TEST ####" | tee std_out.txt
#
#
cd ..
#
#
mv test_inverse/ outputs/
#
#
exit 0
#
# END OF SCRIPT

