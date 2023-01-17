#!/bin/bash
#
mpirun -np 5 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=mt4_valgrind_fwd_mpi5.txt ./ModEM -f -m ../inputs/JMult_test/pr.ws -d ../inputs/JMult_test/de.txt -c ../docs/control_file_template -pd OO_MT4_PREDDATA_MPI5.dat
#
mpirun -np 5 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=mt4_valgrind_j_mpi5.txt ./ModEM -j -m ../inputs/JMult_test/pr.ws -pm ../inputs/JMult_test/dm.ws -d ../inputs/JMult_test/de.txt -c ../docs/control_file_template -pd OO_MT4_JHMAT_MPI5.dat
#
mpirun -np 5 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=mt4_valgrind_jt_mpi5.txt ./ModEM -jt -m ../inputs/JMult_test/pr.ws -d ../inputs/JMult_test/de.txt -c ../docs/control_file_template -dm OO_MT4_DSIGMA_MPI5.rho
#