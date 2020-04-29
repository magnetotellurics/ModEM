#!/bin/bash
#
# ARGUMENTS: 1 - Mod3DMT PATH, 2 - CTRL FILE, 3 - NUMER OF CORES
#
EXEC	=$1
CTRL	=$2
ncores	=$3
#
NAME="${EXEC##*/}"
#
# CREATE TEST OUTPUT FOLDER
mkdir -p test_forward_modelling_$NAME
#
# ENTER TEST OUTPUT FOLDER
cd test_forward_modelling_$NAME/
#
#
echo "#### START FORWARD MODELLING $NAME MPI TEST WITH $ncores CORES AT $now ####" | tee std_out.txt
#
#
echo "#### COMMAND LINE: [mpirun -n $ncores ../$EXEC -W ../$CTRL -v full]" | tee -a std_out.txt
#
#
mpirun -n $ncores ../$EXEC -W ../$CTRL -v full | tee -a std_out.txt
#
# CATCH RESULT
result=$?
#
# TEST RESULT
if [ "$result" -ne "0" ]; then
	#
	#
	echo "TEST FORWARD MODELLING $NAME FAIL: $result" | tee -a std_out.txt
	#
	#
	cd ..
	#
	#
	exit $result
fi
#
#
echo "#### FINISH FORWARD MODELLING $EXEC MPI TEST ####" | tee -a std_out.txt
#
#
cd ..
#
mv test_forward_modelling_$NAME/ outputs/temp/
#
#
exit 0
#
# END OF SCRIPT

