#!/bin/bash
#
mpirun -np 5 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=mt4_valgrind_fwd_mpi5.txt ./ModEM -f -m ../inputs/JMult_test/pr.ws -d ../inputs/JMult_test/de.txt -c ../docs/control_file_template -pd OO_MT4_PREDDATA_MPI5.dat
#
mpirun -np 5 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=mt4_valgrind_j_mpi5.txt ./ModEM -j -m ../inputs/JMult_test/pr.ws -pm ../inputs/JMult_test/dm.ws -d ../inputs/JMult_test/de.txt -c ../docs/control_file_template -pd OO_MT4_JHMAT_MPI5.dat
#
mpirun -np 5 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=mt4_valgrind_jt_mpi5.txt ./ModEM -jt -m ../inputs/JMult_test/pr.ws -d ../inputs/JMult_test/de.txt -c ../docs/control_file_template -dm OO_MT4_DSIGMA_MPI5.rho
#
mpirun -np 3 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=mt4_valgrind_fwd_mpi3.txt ./ModEM -f -m ../inputs/JMult_test/pr.ws -d ../inputs/JMult_test/de.txt -c ../docs/control_file_template -pd OO_MT4_PREDDATA_MPI3.dat
#
mpirun -np 3 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=mt4_valgrind_j_mpi3.txt ./ModEM -j -m ../inputs/JMult_test/pr.ws -pm ../inputs/JMult_test/dm.ws -d ../inputs/JMult_test/de.txt -c ../docs/control_file_template -pd OO_MT4_JHMAT_MPI3.dat
#
mpirun -np 3 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=mt4_valgrind_jt_mpi3.txt ./ModEM -jt -m ../inputs/JMult_test/pr.ws -d ../inputs/JMult_test/de.txt -c ../docs/control_file_template -dm OO_MT4_DSIGMA_MPI3.rho
#
mpirun -np 2 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=mt4_valgrind_fwd_mpi2.txt ./ModEM -f -m ../inputs/JMult_test/pr.ws -d ../inputs/JMult_test/de.txt -c ../docs/control_file_template -pd OO_MT4_PREDDATA_MPI2.dat
#
mpirun -np 2 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=mt4_valgrind_j_mpi2.txt ./ModEM -j -m ../inputs/JMult_test/pr.ws -pm ../inputs/JMult_test/dm.ws -d ../inputs/JMult_test/de.txt -c ../docs/control_file_template -pd OO_MT4_JHMAT_MPI2.dat
#
mpirun -np 2 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=mt4_valgrind_jt_mpi2.txt ./ModEM -jt -m ../inputs/JMult_test/pr.ws -d ../inputs/JMult_test/de.txt -c ../docs/control_file_template -dm OO_MT4_DSIGMA_MPI2.rho
#