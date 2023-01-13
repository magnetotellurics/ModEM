#!/bin/bash
#
valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=mt1_valgrind_jt_serial.txt ./TestSerial -jt -m ../inputs/JMult_test/pr.ws -d ../inputs/JMult_test/d_1tx_1rx.txt -dm OO_DSIGMA.rho
#
valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=mt1_valgrind_j_serial.txt ./TestSerial -j -m ../inputs/JMult_test/pr.ws -pm ../inputs/JMult_test/dm.ws -d ../inputs/JMult_test/d_1tx_1rx.txt -gd JmHat.txt
