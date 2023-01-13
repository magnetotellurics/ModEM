#!/bin/bash
#
# MT nTX=1, nRx=198
#
./TestSerial -f -m ../inputs/JMult_test/pr.ws -d ../inputs/JMult_test/de.txt -pd OO_FWD_PRED_DATA_MT4.dat -c ../docs/control_file_template
#
./TestSerial -jt -m ../inputs/JMult_test/pr.ws -d ../inputs/JMult_test/de.txt -dm OO_JT_DSIGMA_DATA_MT4.rho -c ../docs/control_file_template
#
./TestSerial -j -m ../inputs/JMult_test/pr.ws -pm ../inputs/JMult_test/dm.ws -d ../inputs/JMult_test/de.txt -jm OO_J_JMHAT_DATA_MT4.dat -c ../docs/control_file_template
#
./TestSerial -i -m ../inputs/JMult_test/pr.ws -d ../inputs/JMult_test/de.txt -c ../docs/control_file_template -o OO_DCG_MT4
#