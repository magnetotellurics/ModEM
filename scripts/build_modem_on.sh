#!/bin/bash
#
# NO ARGUMENTS: 
#
# STRING NOW
now=$(date "+%Y/%m/%d - %H:%M:%S")
#
#
cd src/
#
#
echo "#### START BUILD STANDART Mod3DMT AT $now ####"
#
# GRANT PERMISSION TO Configure.3D_MT.OSU.GFortran
chmod 777 CONFIG/Configure.3D_MT.OSU.GFortran
#
# CREATE Makefile_STD
./CONFIG/Configure.3D_MT.OSU.GFortran Makefile_STD MPI
#
# REMOVE OLD OBJECTS
rm -rf ../../objs/
#
# BUILD WITH PLAIN MAKE
make -f Makefile_STD
#
# RENAME EXE
mv Mod3DMT ../bin/Mod3DMT_STD
#
# CATCH RESULT
result=$?
#
# TEST RESULT
if [ "$result" -ne "0" ]; then
	#
	#
	echo "BUILD Mod3DMT_STD FAIL: $result"
	#
	#
	cd ..
	#
	#
	exit $result
fi
#
#
echo "#### FINISH BUILD Mod3DMT_STD QMR ####"
#
#
cd ..
#
#
exit 0
#
#
# END OF SCRIPT

