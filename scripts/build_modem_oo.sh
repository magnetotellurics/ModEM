#!/bin/bash
#
# NO ARGUMENTS: 
#
# CREATE BIN FOLDER (IF NOT EXISTS)
mkdir -p bin
#
# CREATE MAIN OUTPUT FOLDER (IF NOT EXISTS)
mkdir -p outputs/temp
#
echo "#### START BUILD ModEM-OO AT $(date "+%Y/%m/%d - %H:%M:%S")" | tee -a outputs/temp/summary.txt
echo "   #" | tee -a outputs/temp/summary.txt
#
T_START=$(date +%s%3N)
#
# REMOVE OLD OBJECTS
rm -rf ../objs
#
#
cd src
   #
   make clean
   #
   # BUILD MAKEFILE
   ./Configure.3D_MT.OSU.GFortran Makefile OO ModEM-OO.f90
   #
   # BUILD WITH PLAIN MAKE
   make
   #
   # CATCH RESULT
   RESULT=$?
   #
   # TEST RESULT
   if [ "$RESULT" -ne "0" ]; then
      #
      echo "   > BUILD ModEM-OO FAIL: $RESULT" | tee -a ../outputs/temp/summary.txt
      T_END=$(date +%s%3N)
      echo "   > Time Spent: $(( ( $T_END - $T_START ) / 1000 )) seconds" | tee -a ../outputs/temp/summary.txt
      echo "   #" | tee -a ../outputs/temp/summary.txt
      #
      cd ..
      #
      exit $RESULT
      #
   fi
   #
   echo "   > BUILD ModEM-OO SUCCESSFULLY: $RESULT" | tee -a ../outputs/temp/summary.txt
   T_END=$(date +%s%3N)
   echo "   > Time Spent: $(( ( $T_END - $T_START ) / 1000 )) seconds" | tee -a ../outputs/temp/summary.txt
   echo "   #" | tee -a ../outputs/temp/summary.txt
   #
   # RENAME EXE
   mv ModEM-OO ../bin
   #
   cd ..
#
# REMOVE NEW OBJECTS
rm -rf ../objs
#
echo "#### FINISH BUILD ModEM-OO AT $(date "+%Y/%m/%d - %H:%M:%S")" | tee -a outputs/temp/summary.txt
#
# RENAME outputs/temp
mv outputs/temp outputs/$(hostname)_$(date "+%Y%m%d_%H%M%S")
#
exit 0
#
# END OF SCRIPT

