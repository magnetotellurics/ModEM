#!/bin/bash
#
# MT nTX=16, nRx=143
#
./TestSerial -f -m ../inputs/1st_Example/rFile_Model -d ../inputs/1st_Example/rFile_Data_MT -pd OO_FWD_PRED_DATA_MT16.dat -c ../docs/control_file_template
#
./TestSerial -jt -m ../inputs/1st_Example/rFile_Model -d ../inputs/1st_Example/rFile_Data_MT -dm OO_JT_DSIGMA_DATA_MT16.rho -c ../docs/control_file_template
#
./TestSerial -j -m ../inputs/1st_Example/rFile_Model -pm OO_JT_DSIGMA_DATA_MT16.rho -d ../inputs/1st_Example/rFile_Data_MT -jm OO_J_JMHAT_DATA_MT16.dat -c ../docs/control_file_template
#
./TestSerial -i -m ../inputs/1st_Example/rFile_Model -d ../inputs/1st_Example/rFile_Data_MT -c ../docs/control_file_template
#