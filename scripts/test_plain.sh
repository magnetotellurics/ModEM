#!/bin/bash
#
# NO ARGUMENTS: 
#
# STRING NOW
now=$(date "+%Y/%m/%d - %H:%M:%S")
#
# GET ENVIROMENT NUMBER OF CORES 
ncores=$(nproc)
#
#
echo "### START Mod3DMT: $ncores CORES, WITH NO ARGUMENTS ###"
#
#
mpirun -n $ncores ./src/Mod3DMT
#
#
echo "### FINISH Mod3DMT ###"
#
# END OF SCRIPT

