#!/bin/bash
#
# MT nTX=16, nRx=143
#
mpirun -np 5 ./Mod3DMT_STD -F ../../modem-oo/inputs/1st_Example/rFile_Model ../../modem-oo/inputs/1st_Example/rFile_Data_MT ON_FWD_PRED_DATA_MT16.dat ON_FWD_SOLN_MT16.bin FWD_para.dat
#
mpirun -np 5 ./Mod3DMT_STD -T ../../modem-oo/inputs/1st_Example/rFile_Model ../../modem-oo/inputs/1st_Example/rFile_Data_MT ON_JT_DSIGMA_MT16.rho FWD_para.dat
#
mpirun -np 5 ./Mod3DMT_STD -M ../../modem-oo/inputs/1st_Example/rFile_Model ON_JT_DSIGMA_MT16.rho ../../modem-oo/inputs/1st_Example/rFile_Data_MT ON_J_JMHAT_MT16.dat FWD_para.dat
#
mpirun -np 5 ./Mod3DMT_STD -I DCG ../../modem-oo/inputs/1st_Example/rFile_Model ../../modem-oo/inputs/1st_Example/rFile_Data_MT ../../modem-oo/inputs/Others/modem_on_inv_control_file.txt FWD_para.dat
#
# CSEM nTX=45, nRx=65
#
mpirun -np 6 ./Mod3DMT_STD -F ../../modem-oo/inputs/Naser_CSEM/rFile_Model ../../modem-oo/inputs/Naser_CSEM/rFile_Data_E1 ON_FWD_PRED_DATA_CSEM45.dat ON_FWD_SOLN_CSEM45.bin FWD_para.dat
#
mpirun -np 6 ./Mod3DMT_STD -T ../../modem-oo/inputs/Naser_CSEM/rFile_Model ../../modem-oo/inputs/Naser_CSEM/rFile_Data_E1 ON_JT_DSIGMA_CSEM45.rho FWD_para.dat
#
mpirun -np 6 ./Mod3DMT_STD -M ../../modem-oo/inputs/Naser_CSEM/rFile_Model ON_JT_DSIGMA_CSEM45.rho ../../modem-oo/inputs/Naser_CSEM/rFile_Data_E1 ON_J_JMHAT_CSEM45.dat FWD_para.dat
#
mpirun -np 6 ./Mod3DMT_STD -I DCG ../../modem-oo/inputs/Naser_CSEM/rFile_Model ../../modem-oo/inputs/Naser_CSEM/rFile_Data_E1 ../../modem-oo/inputs/Others/modem_on_inv_control_file.txt FWD_para.dat
#