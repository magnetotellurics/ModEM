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
echo "### START BUILD Mod3DMT QMR AT $now ###" >> build_modem.txt
#
# DEFINE THE SELA SOLVING METHOD: QMR or BICG
chmod 777 Configure.3D_MT.GFortran &>> build_modem.txt
#
#
./Configure.3D_MT.GFortran Makefile MPI &>> build_modem.txt
#
# REMOVE OLD OBJECTS
rm -rf objs/3D_MT/GFortReleaseMPI/* &>> build_modem.txt
#
# BUILD WITH PLAIN MAKE
make &>> build_modem.txt
#
#
echo "### FINISH BUILD Mod3DMT QMR ###" >> build_modem.txt
#
#
mv build_modem.txt ../outputs/
#
#
cd ..
#
#
# END OF SCRIPT

