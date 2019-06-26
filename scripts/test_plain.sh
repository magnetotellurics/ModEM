#!/bin/bash
#
# NO ARGUMENTS: 
#
# STRING NOW
now=$(date "+%Y_%m_%d_%H_%M_%S")
#
# GET ENVIROMENT NUMBER OF CORES 
ncores=$(nproc)
#
#
echo "### START Mod3DMT: $ncores CORES, WITH NO ARGUMENTS ###" >> test_plain_$now.txt
#
#
mpirun -n $ncores ./src/Mod3DMT &>> test_plain_$now.txt
#
#
echo "### FINISH Mod3DMT ###" >> test_plain_$now.txt
#
#
mv test_plain_$now.txt outputs/
#
#
# END OF SCRIPT

