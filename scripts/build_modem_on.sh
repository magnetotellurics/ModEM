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
########################
# STANDART VERSION
########################
#
echo "#### START BUILD STANDART Mod3DMT AT $now ####"
#
T_START=$(date +%s%3N)
#
echo "	> START BUILD MODEM STD" | tee -a ../outputs/temp/summary.txt
#
# GRANT PERMISSION TO Configure.3D_MT.OSU.GFortran
chmod 777 CONFIG/Configure.3D_MT_EM1DCSEM.OSU.GFortran
#
# CREATE Makefile_STD
./CONFIG/Configure.3D_MT_EM1DCSEM.OSU.GFortran Makefile_STD MPI
#
# REMOVE OLD OBJECTS
rm -rf ../../objs/
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
mv Mod3DMT ../bin/ModEM_STD
#
########################
# SP2 VERSION
########################
#
echo "#### START BUILD SP2 Mod3DMT AT $now ####"
#
T_START=$(date +%s%3N)
#
echo "	> START BUILD MODEM SP2" | tee -a ../outputs/temp/summary.txt
#
# GRANT PERMISSION TO Configure.3D_MT_SP2.OSU.GFortran
chmod 777 CONFIG/Configure.3D_MT_SP2_EM1DCSEM.OSU.GFortran
#
# CREATE Makefile_STD
./CONFIG/Configure.3D_MT_SP2_EM1DCSEM.OSU.GFortran Makefile_SP2 MPI
#
# REMOVE OLD OBJECTS
rm -rf ../../objs/
#
# BUILD WITH PLAIN MAKE
make -f Makefile_SP2
#
# CATCH RESULT
RESULT=$?
#
# TEST RESULT
if [ "$RESULT" -ne "0" ]; then
	#
	echo "	> BUILD MODEM SP2 FAIL: $RESULT" | tee -a ../outputs/temp/summary.txt
	T_END=$(date +%s%3N)
	echo "	> Time Spent: $(( ( $T_END - $T_START ) / 1000 )) seconds" | tee -a ../outputs/temp/summary.txt
	echo "	#" | tee -a ../outputs/temp/summary.txt
	#
	cd ..
	#
	exit $RESULT
fi
#
echo "#### FINISH BUILD Mod3DMT SP2 ####"
#
echo "	> BUILD MODEM SP2 PASS: $RESULT" | tee -a ../outputs/temp/summary.txt
T_END=$(date +%s%3N)
echo "	> Time Spent: $(( ( $T_END - $T_START ) / 1000 )) seconds" | tee -a ../outputs/temp/summary.txt
echo "	#" | tee -a ../outputs/temp/summary.txt
#
# RENAME EXE
mv Mod3DMT ../bin/ModEM_SP2
#
# REMOVE NEW OBJECTS
rm -rf objs/
cd ..
#
exit 0
#
# END OF SCRIPT

