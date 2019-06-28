#!/bin/bash
#
# ARGUMENTS: 1 - Mod3DMT EXECUTABLE, 2 - MODEL FILE, 3 - DATA FILE, 4 - dMODEL FILE
EXEC=$1
MODEL=$2
DATA=$3
dMODEL=$4
#
# STRING NOW
now=$(date "+%Y/%m/%d - %H:%M:%S")
#
# GET ENVIROMENT NUMBER OF CORES 
ncores=$(nproc)
#
# REMOVE TEST OUTPUT FOLDER FROM MAIN OUTPUT FOLDER
rm -rf outputs/test_grad/
#
# CREATE TEST OUTPUT FOLDER
mkdir -p test_grad
#
# ENTER TEST OUTPUT FOLDER
cd test_grad/
#
#
echo "#### START GRAD MPI TEST WITH $ncores CORES AT $now ####" >> std_out.txt
#
#
echo "#### COMMAND LINE: [mpirun -n $ncores ../$EXEC -g ../$MODEL ../$DATA ../$dMODEL -v full]" >> std_out.txt
#
#
mpirun -n $ncores ../$EXEC -g ../$MODEL ../$DATA ../$dMODEL -v full &>> std_out.txt
#
# CATCH RESULT
result=$?
#
# TEST RESULT
if [ "$result" -ne "0" ]; then
	#
	#
	echo "TEST GRAD FAIL: $result"
	#
	#
	cd ..
	#
	#
	exit $result
fi
#
#
echo "#### FINISH GRAD MPI TEST ####" >> std_out.txt
#
#
cd ..
#
#
mv test_grad/ outputs/
#
#
exit 0
#
# END OF SCRIPT

