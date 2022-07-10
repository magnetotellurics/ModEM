#!/bin/bash
#
#--suppressions=../../../MEETING_INPUTS/modem_oo_suppresions.supp
#
mpirun -np 2 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=oo_mt_valgrind_fwd_mpi.txt ../test_FWD_MPI --forward --model ../../../MEETING_INPUTS/rFile_Model  --data ../../../MEETING_INPUTS/rFile_Data_MT --control ../../inputs/Others/first_control_file.txt | tee -a mt_output.txt
#
mpirun -np 2 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=oo_csem_valgrind_fwd_mpi.txt ../test_FWD_MPI --forward --model ../../../MEETING_INPUTS/rFile_Model  --data ../../../MEETING_INPUTS/rFile_Data_CSEM --control ../../inputs/Others/first_control_file.txt | tee -a output.txt
#
