#!/bin/bash
#
# ARGUMENTS: 1 - Mod3DMT EXECUTABLE, 2 - CTRL FILE, 3 - NUMER OF CORES
EXEC=$1
CTRL=$2
ncores=$3
#
# STRING NOW
now=$(date "+%Y/%m/%d - %H:%M:%S")
#
exec_name="${EXEC##*/}"
#
# CREATE TEST OUTPUT FOLDER
mkdir -p test_forward_modelling_$exec_name
#
# ENTER TEST OUTPUT FOLDER
cd test_forward_modelling_$exec_name/
#
#
echo "#### START FORWARD $exec_name MPI TEST WITH $ncores CORES AT $now ####" | tee std_out.txt
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
	echo "TEST FORWARD $exec_name FAIL: $result" | tee -a std_out.txt
	#
	#
	cd ..
	#
	#
	exit $result
fi
#
#
echo "#### FINISH FORWARD $EXEC MPI TEST ####" | tee -a std_out.txt
#
#
cd ..
#
# BUILD bin/SolverDiagnostic3D
cd tools/SolverDiagnostic3D/
bash build_linux.sh
#
cd ../../test_forward_modelling_$exec_name
../tools/SolverDiagnostic3D/bin/SolverDiagnostic3D QMR* ../tools/MathBox/mathbox-bundle.js
cd ..
#
#
mv test_forward_modelling_$exec_name/ outputs/temp/
#
#
exit 0
#
# END OF SCRIPT

