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
echo "#### START BUILD Mod3DMT AT $now ####"
#
# 
chmod 777 CONFIG/Configure.3D_MT.OSU.GFortran
#
#
./CONFIG/Configure.3D_MT.OSU.GFortran Makefile MPI
#
# REMOVE OLD OBJECTS
rm -rf objs/3D_MT/GFortReleaseMPI/*
#
# BUILD WITH PLAIN MAKE
make
#
# CATCH RESULT
result=$?
#
# TEST RESULT
if [ "$result" -ne "0" ]; then
	#
	#
	echo "build_modem_on FAIL: $result"
	#
	#
	cd ..
	#
	#
	exit $result
fi
#
#
echo "#### FINISH BUILD Mod3DMT QMR ####"
#
#
cd ..
#
#
exit 0
#
#
# END OF SCRIPT

