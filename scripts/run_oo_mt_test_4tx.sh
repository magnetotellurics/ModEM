#!/bin/bash
#
# SERIAL MT4
#
# FWD: 5.734s, 8.297s
#
./ModEM_SERIAL -f -m ../inputs/JMult_test/pr.ws -d ../inputs/JMult_test/de.txt -c ../docs/control_file_template -pd OO_PREDDATA_MT4_SERIAL.dat
#
# JMult_t: 23.578s, 39.578s
#
./ModEM_SERIAL -jt -m ../inputs/JMult_test/pr.ws -d ../inputs/JMult_test/de.txt -c ../docs/control_file_template -dm OO_DSIGMA_MT4_SERIAL.rho
#
# JMult: 23.500s, 37.531s
#
./ModEM_SERIAL -j -m ../inputs/JMult_test/pr.ws -pm ../inputs/JMult_test/dm.ws -d ../inputs/JMult_test/de.txt -c ../docs/control_file_template -jm OO_JMHAT_SERIAL.dat
#
# DCG: 
#
./ModEM_SERIAL -i -m ../inputs/JMult_test/pr.ws -d ../inputs/JMult_test/de.txt -c ../docs/control_file_template -o OO_DCG_SERIAL
#
# PARALLEL MT4
#
# FWD: 6.750s
#
mpirun -np 5 ./ModEM_MPI -f -m ../inputs/JMult_test/pr.ws -d ../inputs/JMult_test/de.txt -c ../docs/control_file_template -pd OO_PREDDATA_MT4_NP5.dat
#
# JMult_t: 14.094s
#
mpirun -np 5 ./ModEM_MPI -jt -m ../inputs/JMult_test/pr.ws -d ../inputs/JMult_test/de.txt -c ../docs/control_file_template -dm OO_DSIGMA_MT4_NP5.rho
#
# JMult: 11.219s
#
mpirun -np 5 ./ModEM_MPI -j -m ../inputs/JMult_test/pr.ws -pm ../inputs/JMult_test/dm.ws -d ../inputs/JMult_test/de.txt -c ../docs/control_file_template -jm OO_JMHAT_MT4_NP5.dat
#
# DCG: 1195.125s 1255.109s
#
mpirun -np 5 ./ModEM_MPI -i -m ../inputs/JMult_test/pr.ws -d ../inputs/JMult_test/de.txt -c ../docs/control_file_template -o DCG_MPI_MT4
#