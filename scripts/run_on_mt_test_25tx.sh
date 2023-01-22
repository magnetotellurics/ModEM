#!/bin/bash
#
# PARALLEL MT: 25 TX, 65 RX
#
# FWD: 
#
mpirun -np 26 --hostfile ../../modem-oo/inputs/Naser_MT/hostfile25_6_nodes.txt ./Mod3DMT_STD -F ../../modem-oo/inputs/Naser_MT/rFile_Model ../../modem-oo/inputs/Naser_MT/rFile_Data_fix ON_FWD_PRED_DATA_MT25.dat ON_FWD_SOLN_MT25.bin control.fwd
#
# JMult_t: 
#
mpirun -np 26 --hostfile ../../modem-oo/inputs/Naser_MT/hostfile25_6_nodes.txt ./Mod3DMT_STD -T ../../modem-oo/inputs/Naser_MT/rFile_Model ../../modem-oo/inputs/Naser_MT/rFile_Data_fix ON_JT_DSIGMA_MT25.rho control.fwd
#
# JMult: 
#
mpirun -np 26 --hostfile ../../modem-oo/inputs/Naser_MT/hostfile25_6_nodes.txt ./Mod3DMT_STD -M ../../modem-oo/inputs/Naser_MT/rFile_Model ../../modem-oo/inputs/Naser_MT/rFile_Model ../../modem-oo/inputs/Naser_MT/rFile_Data_fix ON_J_JMHAT_MT25.dat control.fwd
#
# DCG: 
#
mpirun -np 26 --hostfile ../../modem-oo/inputs/Naser_MT/hostfile25_6_nodes.txt ./Mod3DMT_STD -I DCG ../../modem-oo/inputs/Naser_MT/rFile_Model ../../modem-oo/inputs/Naser_MT/rFile_Data_fix control.inv control.fwd
#