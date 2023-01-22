#!/bin/bash
#
# PARALLEL MT: 25 TX, 65 RX
#
# FWD: 
#
mpirun -np 26 --hostfile ../../hostfile_mt ./ModEM_MPI -f -m ../inputs/Naser_MT/rFile_Model -d ../inputs/Naser_MT/rFile_Data_fix -c ../docs/control_file_template -pd OO_PREDDATA_MT25_MPI.dat
#
# JMult_t: 
#
mpirun -np 26 --hostfile ../../hostfile_mt ./ModEM_MPI -jt -m ../inputs/Naser_MT/rFile_Model -d ../inputs/Naser_MT/rFile_Data_fix -c ../docs/control_file_template -dm OO_DSIGMA_MT25_MPI.rho
#
# JMult: 
#
mpirun -np 26 --hostfile ../../hostfile_mt ./ModEM_MPI -j -m ../inputs/Naser_MT/rFile_Model -pm ../inputs/Naser_MT/rFile_Model -d ../inputs/Naser_MT/rFile_Data_fix -c ../docs/control_file_template -jm OO_JMHAT_MT25_MPI.dat
#
# DCG: 
#
mpirun -np 26 --hostfile ../../hostfile_mt ./ModEM_MPI -i -m ../inputs/Naser_MT/rFile_Model -d ../inputs/Naser_MT/rFile_Data_fix -c ../docs/control_file_template -o DCG_MT25_MPI
#
