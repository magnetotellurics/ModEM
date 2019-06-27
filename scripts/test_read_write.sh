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
rm -rf outputs/test_read_write/
#
# CREATE TEST OUTPUT FOLDER
mkdir -p test_read_write
#
# ENTER TEST OUTPUT FOLDER
cd test_read_write/
#
#
echo "#### START READ_WRITE MPI TEST WITH $ncores CORES AT $now ####\n" >> std_out.txt
#
#
echo "#### COMMAND LINE: [mpirun -n $ncores ../$EXEC -R ../$MODEL ../$DATA wFile_Model.sw wFile_Data.dat -v full]" >> std_out.txt
#
#
mpirun -n $ncores ../$EXEC -R ../$MODEL ../$DATA wFile_Model.sw wFile_Data.dat -v full &>> std_out.txt
#
# CATCH RESULT
result=$?
#
# TEST RESULT
if [ "$result" -ne "0" ]; then
	#
	#
	echo "TEST READ_WRITE FAIL: $result"
	#
	#
	cd ..
	#
	#
	exit $result
fi
#
#
echo "### FINISH READ_WRITE MPI TEST ###\n" >> std_out.txt
#
#
cd ..
#
#
mv test_read_write/ outputs/
#
#
exit 0
#
# END OF SCRIPT

