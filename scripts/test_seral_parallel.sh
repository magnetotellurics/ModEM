#!/bin/bash
#
# SERIAL MT4
#
# FWD: 5.734s
#
./ModEM_SERIAL -f -m ../inputs/JMult_test/pr.ws -d ../inputs/JMult_test/de.txt -c ../docs/control_file_template -pd OO_PREDDATA_MT4_SERIAL.dat
#
# JMult_t: 23.578s
#
./ModEM_SERIAL -jt -m ../inputs/JMult_test/pr.ws -d ../inputs/JMult_test/de.txt -c ../docs/control_file_template -dm OO_DSIGMA_MT4_SERIAL.rho
#
# JMult: 23.500s
#
./ModEM_SERIAL -j -m ../inputs/JMult_test/pr.ws -pm ../inputs/JMult_test/dm.ws -d ../inputs/JMult_test/de.txt -c ../docs/control_file_template -jm OO_JMHAT_SERIAL.dat
#
# PARALLEL MT4
#
# FWD: 6.750s
#
mpirun -np 5 ./ModEM_MPI -f -m ../inputs/JMult_test/pr.ws -d ../inputs/JMult_test/de.txt -c ../docs/control_file_template -pd OO_PREDDATA_MT4_NP5.dat
#
# JMult_t: 14.094s
#
mpirun -np 5 ./ModEM_MPI -jt -m ../inputs/JMult_test/pr.ws -d ../inputs/JMult_test/de.txt -c ../docs/control_file_template -dm OO_DSIGMA_MT4.rho
#
# JMult: 11.219s
#
mpirun -np 5 ./ModEM_MPI -j -m ../inputs/JMult_test/pr.ws -pm ../inputs/JMult_test/pr.ws -d ../inputs/JMult_test/de.txt -c ../docs/control_file_template -jm OO_JMHAT_MT4_NP5.dat
#