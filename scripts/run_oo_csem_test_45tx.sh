#!/bin/bash
#
# PARALLEL CSEM: 45 TX, 65 RX
#
# FWD: 
#
mpirun -np 46 --hostfile ../inputs/Naser_CSEM/hostfile46_10_nodes.txt ./ModEM_MPI -f -m ../inputs/Naser_CSEM/rFile_Model -d ../inputs/Naser_CSEM/rFile_Data_fix -c ../docs/control_file_template -pd OO_PREDDATA_CSEM45_MPI.dat
#
# JMult_t: 
#
mpirun -np 46 --hostfile ../inputs/Naser_CSEM/hostfile46_10_nodes.txt ./ModEM_MPI -jt -m ../inputs/Naser_CSEM/rFile_Model -d ../inputs/Naser_CSEM/rFile_Data_fix -c ../docs/control_file_template -dm OO_DSIGMA_CSEM45_MPI.rho
#
# JMult: 
#
mpirun -np 46 --hostfile ../inputs/Naser_CSEM/hostfile46_10_nodes.txt ./ModEM_MPI -j -m ../inputs/Naser_CSEM/rFile_Model -pm ../inputs/Naser_CSEM/rFile_Model -d ../inputs/Naser_CSEM/rFile_Data_fix -c ../docs/control_file_template -jm OO_JMHAT_CSEM45_MPI.dat
#
# DCG: 
#
mpirun -np 46 --hostfile ../inputs/Naser_CSEM/hostfile46_10_nodes.txt ./ModEM_MPI -i -m ../inputs/Naser_CSEM/rFile_Model -d ../inputs/Naser_CSEM/rFile_Data_fix -c ../docs/control_file_template -o DCG_CSEM45_MPI
#