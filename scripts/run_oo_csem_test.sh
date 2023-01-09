#!/bin/bash
#
# CSEM nTX=45, nRx=65
#
./TestSerial -f -m ../inputs/Naser_CSEM/rFile_Model -d ../inputs/Naser_CSEM/rFile_Data_E1 -pd OO_FWD_PRED_DATA_CSEM45.dat -c ../docs/control_file_template
#
./TestSerial -jt -m ../inputs/Naser_CSEM/rFile_Model -d ../inputs/Naser_CSEM/rFile_Data_E1 -dm OO_JT_DSIGMA_DATA_CSEM45.rho -c ../docs/control_file_template
#
./TestSerial -j -m ../inputs/Naser_CSEM/rFile_Model -pm OO_JT_DSIGMA_DATA_CSEM45.rho -d ../inputs/Naser_CSEM/rFile_Data_E1 -jm OO_J_JMHAT_DATA_CSEM45.dat -c ../docs/control_file_template
#
./TestSerial -i -m ../inputs/Naser_CSEM/rFile_Model -d ../inputs/Naser_CSEM/rFile_Data_E1 -c ../docs/control_file_template
#