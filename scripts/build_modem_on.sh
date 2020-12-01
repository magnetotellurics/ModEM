#!/bin/bash
#
# NO ARGUMENTS: 
mkdir -p bin
#
# STRING NOW
now=$(date "+%Y/%m/%d - %H:%M:%S")
#
cd src/
#
# REMOVE OLD OBJECTS
rm -rf objs/
#
#
echo "#### START BUILD STANDART Mod3DMT AT $now ####"
#
T_START=$(date +%s%3N)
#
echo "	> START BUILD MODEM STD" | tee -a ../outputs/temp/summary.txt
#
# GRANT PERMISSION TO Configure.3D_MT.OSU.GFortran
chmod 777 CONFIG/Configure.3D_MT.OSU.GFortran
#
# CREATE Makefile_STD
./CONFIG/Configure.3D_MT.OSU.GFortran Makefile_STD MPI
#
# BUILD WITH PLAIN MAKE
make -f Makefile_STD
#
# CATCH RESULT
RESULT=$?
#
# TEST RESULT
if [ "$RESULT" -ne "0" ]; then
	#
	echo "	> BUILD MODEM STD FAIL: $RESULT" | tee -a ../outputs/temp/summary.txt
	T_END=$(date +%s%3N)
	echo "	> Time Spent: $(( ( $T_END - $T_START ) / 1000 )) seconds" | tee -a ../outputs/temp/summary.txt
	echo "	#" | tee -a ../outputs/temp/summary.txt
	#
	cd ..
	#
	exit $RESULT
fi
#
echo "#### FINISH BUILD Mod3DMT STD ####"
#
echo "	> BUILD MODEM STD PASS: $RESULT" | tee -a ../outputs/temp/summary.txt
T_END=$(date +%s%3N)
echo "	> Time Spent: $(( ( $T_END - $T_START ) / 1000 )) seconds" | tee -a ../outputs/temp/summary.txt
echo "	#" | tee -a ../outputs/temp/summary.txt
#
# RENAME EXE
mv Mod3DMT ../bin/ModEM_baseline
#
# REMOVE NEW OBJECTS
rm -rf objs/
#
#
cd ..
#
exit 0
#
# END OF SCRIPT

