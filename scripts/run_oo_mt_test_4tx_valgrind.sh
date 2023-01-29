#!/bin/bash
#
# SERIAL MT4
#
# FWD: 
#
valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=mt4_valgrind_fwd_serial.txt ./ModEM_SERIAL -f -m ../inputs/JMult_test/pr.ws -d ../inputs/JMult_test/de.txt -c ../docs/control_file_template -pd OO_PREDDATA_MT4_SERIAL.dat
#
# JMult_t: 
#
valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=mt4_valgrind_jt_serial.txt ./ModEM_SERIAL -jt -m ../inputs/JMult_test/pr.ws -d ../inputs/JMult_test/de.txt -c ../docs/control_file_template -dm OO_DSIGMA_MT4_SERIAL.rho
#
# JMult: 
#
valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=mt4_valgrind_j_serial.txt ./ModEM_SERIAL -j -m ../inputs/JMult_test/pr.ws -pm ../inputs/JMult_test/dm.ws -d ../inputs/JMult_test/de.txt -c ../docs/control_file_template -jm OO_JMHAT_SERIAL.dat
#
# DCG: 
#
valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=mt4_valgrind_dcg_serial.txt ./ModEM_SERIAL -i -m ../inputs/JMult_test/pr.ws -d ../inputs/JMult_test/de.txt -c ../docs/control_file_template -o OO_DCG_SERIAL
#
# PARALLEL MT4
#
# FWD: 
#
mpirun -np 5 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=mt4_valgrind_fwd_mpi5.txt ./ModEM_MPI -f -m ../inputs/JMult_test/pr.ws -d ../inputs/JMult_test/de.txt -c ../docs/control_file_template -pd OO_PREDDATA_MT4_NP5.dat
#
# JMult_t: 
#
mpirun -np 5 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=mt4_valgrind_jt_mpi5.txt ./ModEM_MPI -jt -m ../inputs/JMult_test/pr.ws -d ../inputs/JMult_test/de.txt -c ../docs/control_file_template -dm OO_DSIGMA_MT4_NP5.rho
#
# JMult: 
#
mpirun -np 5 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=mt4_valgrind_j_mpi5.txt ./ModEM_MPI -j -m ../inputs/JMult_test/pr.ws -pm ../inputs/JMult_test/dm.ws -d ../inputs/JMult_test/de.txt -c ../docs/control_file_template -jm OO_JMHAT_MT4_NP5.dat
#
# DCG: 
#
mpirun -np 5 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=mt4_valgrind_dcg_mpi5.txt ./ModEM_MPI -i -m ../inputs/JMult_test/pr.ws -d ../inputs/JMult_test/de.txt -c ../docs/control_file_template -o DCG_MPI_MT4
#