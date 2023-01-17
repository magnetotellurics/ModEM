#!/bin/bash
#
# CSEM nTX=1, nRx=13
#
./TestSerial -f -m ../inputs/1st_Example/rFile_Model -d ../inputs/1st_Example/rFile_Data_CSEM -pd OO_FWD_PRED_DATA_CSEM1.dat -c ../docs/control_file_template
#
./TestSerial -jt -m ../inputs/1st_Example/rFile_Model -d ../inputs/1st_Example/rFile_Data_CSEM -dm OO_JT_DSIGMA_DATA_CSEM1.rho -c ../docs/control_file_template
#
./TestSerial -j -m ../inputs/1st_Example/rFile_Model -pm OO_JT_DSIGMA_DATA_CSEM1.rho -d ../inputs/1st_Example/rFile_Data_CSEM -jm OO_J_JMHAT_DATA_CSEM1.dat -c ../docs/control_file_template
#
#./TestSerial -i -m ../inputs/1st_Example/rFile_Model -d ../inputs/1st_Example/rFile_Data_CSEM -c ../docs/control_file_template
#