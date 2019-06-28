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
rm -rf outputs/test_forward/
#
# CREATE TEST OUTPUT FOLDER
mkdir -p test_forward
#
# ENTER TEST OUTPUT FOLDER
cd test_forward/
#
#
echo "#### START FORWARD MPI TEST WITH $ncores CORES AT $now ####" >> std_out.txt
#
#
echo "#### COMMAND LINE: [mpirun -n $ncores ../$EXEC -F ../$MODEL ../$DATA wFile_Data.dat wFile_EMsoln -v full]" >> std_out.txt
#
#
mpirun -n $ncores ../$EXEC -F ../$MODEL ../$DATA wFile_Data.dat wFile_EMsoln -v full &>> std_out.txt
#
# CATCH RESULT
result=$?
#
# TEST RESULT
if [ "$result" -ne "0" ]; then
	#
	#
	echo "TEST FORWARD FAIL: $result"
	#
	#
	cd ..
	#
	#
	exit $result
fi
#
#
echo "#### FINISH FORWARD MPI TEST ####" >> std_out.txt
#
#
cd ..
#
#
mv test_forward/ outputs/
#
#
exit 0
#
# END OF SCRIPT

