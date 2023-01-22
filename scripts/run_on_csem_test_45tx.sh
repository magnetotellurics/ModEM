#
# PARALLEL CSEM: 45 TX, 65 RX
#
# FWD: 
#
mpirun -np 46 --hostfile ../../modem-oo/inputs/Naser_CSEM/hostfile46_10_nodes.txt ./Mod3DMT_STD -F ../../modem-oo/inputs/Naser_CSEM/rFile_Model ../../modem-oo/inputs/Naser_CSEM/rFile_Data_fix ON_FWD_PRED_DATA_CSEM45.dat ON_FWD_SOLN_CSEM45.bin control.fwd
#
# JMult_t: 
#
mpirun -np 46 --hostfile ../../modem-oo/inputs/Naser_CSEM/hostfile46_10_nodes.txt ./Mod3DMT_STD -T ../../modem-oo/inputs/Naser_CSEM/rFile_Model ../../modem-oo/inputs/Naser_CSEM/rFile_Data_fix ON_JT_DSIGMA_CSEM45.rho control.fwd
#
# JMult: 
#
mpirun -np 46 --hostfile ../../modem-oo/inputs/Naser_CSEM/hostfile46_10_nodes.txt ./Mod3DMT_STD -M ../../modem-oo/inputs/Naser_CSEM/rFile_Model ../../modem-oo/inputs/Naser_CSEM/rFile_Model ../../modem-oo/inputs/Naser_CSEM/rFile_Data_fix ON_J_JMHAT_CSEM45.dat control.fwd
#
# DCG: 
#
mpirun -np 46 --hostfile ../../modem-oo/inputs/Naser_CSEM/hostfile46_10_nodes.txt ./Mod3DMT_STD -I DCG ../../modem-oo/inputs/Naser_CSEM/rFile_Model ../../modem-oo/inputs/Naser_CSEM/rFile_Data_fix control.inv control.fwd
#