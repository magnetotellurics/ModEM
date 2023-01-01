#!/bin/bash
#
mpirun -np 33 ./Mod3DMT_STD -F ../../modem-oo/inputs/1st_Example/rFile_Model_homo10ohm ../../modem-oo/inputs/1st_Example/rFile_Data_MT fwd_mt4_pred_data.txt fwd_mt4_esolution.bin ../../modem-oo/inputs/Others/modem_on_control_file.txt
#
mpirun -np 33 ./Mod3DMT_STD -T ../../modem-oo/inputs/1st_Example/rFile_Model_homo10ohm ../../modem-oo/inputs/1st_Example/rFile_Data_MT adjt_mt4_dsigma ../../modem-oo/inputs/Others/modem_on_control_file.txt
#
mpirun -np 33 ./Mod3DMT_STD -M ../../modem-oo/inputs/1st_Example/rFile_Model_homo10ohm ../../modem-oo/inputs/1st_Example/rFile_Model_trim ../../modem-oo/inputs/1st_Example/rFile_Data_MT mt4_adj_grad_data_mpi.txt ../../modem-oo/inputs/Others/modem_on_control_file.txt
#
mpirun -np 46 ./Mod3DMT_STD -F ../../modem-oo/inputs/Naser_CSEM/rFile_Model ../../modem-oo/inputs/Naser_CSEM/rFile_Data_fix fwd_csem_pred_data.txt fwd_csem_esolution.bin ../../modem-oo/inputs/Others/modem_on_control_file.txt
#
mpirun -np 46 ./Mod3DMT_STD -T ../../modem-oo/inputs/Naser_CSEM/rFile_Model ../../modem-oo/inputs/Naser_CSEM/rFile_Data_fix adjt_csem_dsigma ../../modem-oo/inputs/Others/modem_on_control_file.txt
#
mpirun -np 46 ./Mod3DMT_STD -M ../../modem-oo/inputs/Naser_CSEM/rFile_Model ../../modem-oo/inputs/Naser_CSEM/rFile_Model ../../modem-oo/inputs/Naser_CSEM/rFile_Data_fix csem_adj_grad_data_mpi.txt ../../modem-oo/inputs/Others/modem_on_control_file.txt
#