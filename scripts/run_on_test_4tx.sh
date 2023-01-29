#!/bin/bash
#
# FWD MT4
#
mpirun -np 5 ./Mod3DMT_STD -F ../../modem-oo/inputs/JMult_test/pr.ws ../../modem-oo/inputs/JMult_test/de.txt ON_PRED_DATA_MT4.dat ON_ESOLN_MT4.bin FWD_para.dat
#
# JMULT_T MT4
#
mpirun -np 5 ./Mod3DMT_STD -T ../../modem-oo/inputs/JMult_test/pr.ws ../../modem-oo/inputs/JMult_test/de.txt ON_DSIGMA_MT4.rho FWD_para.dat
#
# JMULT MT4
#
mpirun -np 5 ./Mod3DMT_STD -M ../../modem-oo/inputs/JMult_test/pr.ws ../../modem-oo/inputs/JMult_test/dm.ws ../../modem-oo/inputs/JMult_test/de.txt ON_JMHAT_MT4.dat FWD_para.dat
#
# INVERSION DCG MT4 : 5.93385410 min OR 355 sec
#
#mpirun -np 5 ./Mod3DMT_STD -I DCG ../../modem-oo/inputs/JMult_test/pr.ws ../../modem-oo/inputs/JMult_test/de.txt ../../modem-oo/inputs/Others/modem_on_inv_control_file.txt FWD_para.dat
#
