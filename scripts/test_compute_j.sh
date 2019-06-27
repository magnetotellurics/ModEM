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
rm -rf outputs/test_compute_j/
#
# CREATE TEST OUTPUT FOLDER
mkdir -p test_compute_j
#
# ENTER TEST OUTPUT FOLDER
cd test_compute_j/
#
#
echo "#### START COMPUT_J MPI TEST WITH $ncores CORES AT $now ####\n" >> std_out.txt
#
#
echo "#### COMMAND LINE: [mpirun -n $ncores ../$EXEC -R ../inputs/rFile_Model.sw ../inputs/rFile_Data.dat wFile_Model.sw wFile_Data.dat -v full]" >> std_out.txt
#
#
mpirun -n $ncores ../$EXEC -J ../$MODEL ../$DATA wFile_Sens -v full &>> std_out.txt
#
# CATCH RESULT
result=$?
#
# TEST RESULT
if [ "$result" -ne "0" ]; then
	#
	#
	echo "TEST COMPUT_J FAIL: $result"
	#
	#
	cd ..
	#
	#
	exit $result
fi
#
#
echo "### FINISH COMPUT_J MPI TEST ###\n" >> std_out.txt
#
#
cd ..
#
#
mv test_compute_j/ outputs/
#
#
exit 0
#
# END OF SCRIPT

