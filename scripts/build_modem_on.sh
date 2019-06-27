#!/bin/bash
#
# NO ARGUMENTS: 
#
# STRING NOW
now=$(date "+%Y/%m/%d - %H:%M:%S")
#
# GENERATE MAIN OUTPUT FOLDER (IF NOT EXISTS)
mkdir -p outputs/
#
#
cd src/
#
#
echo "### START BUILD Mod3DMT AT $now ###" #>> build_modem.txt
#
# 
chmod 777 Configure.3D_MT.GFortran #&>> build_modem.txt
#
#
./Configure.3D_MT.GFortran Makefile MPI OUT #&>> build_modem.txt
#
# REMOVE OLD OBJECTS
rm -rf objs/3D_MT/GFortReleaseMPI/* #&>> build_modem.txt
#
# BUILD WITH PLAIN MAKE
make #&>> build_modem.txt
#
# CATCH RESULT
result=$?
#
# TEST RESULT
if [ "$result" -ne "0" ]; then
   echo "build_modem_on FAIL: $result" #&>> build_modem.txt
   exit $result
fi
#
#
echo "### FINISH BUILD Mod3DMT QMR ###" #>> build_modem.txt
#
#
#mv build_modem.txt ../outputs/
#
#
cd ..
#
exit 0
#
#
# END OF SCRIPT

