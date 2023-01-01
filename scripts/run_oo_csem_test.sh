#!/bin/bash
# #
# # MPI VALGRIND HUGE
# #
# mpirun -np 2 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=csem_valgrind_fwd_mpi.txt ./TestMPI -f -m ../inputs/Naser_CSEM/rFile_Model -d ../inputs/Naser_CSEM/rFile_Data_fix -pd pred_data_fwd_csem_mpi.txt -es esolution_fwd_csem_mpi.bin -c ../docs/control_file_template
# #
# mpirun -np 2 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=csem_valgrind_adjt_mpi.txt ./TestMPI -jt -m ../inputs/Naser_CSEM/rFile_Model -d ../inputs/Naser_CSEM/rFile_Data_fix -pd pred_data_adjt_csem_mpi.txt -dm adjt_csem_sigma_mpi.txt -es esolution_adjt_csem_mpi.bin -c ../docs/control_file_template
# #
# mpirun -np 2 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=csem_valgrind_adj_mpi.txt ./TestMPI -j -m ../inputs/Naser_CSEM/rFile_Model -pm adjt_csem_sigma_mpi.txt -d ../inputs/Naser_CSEM/rFile_Data_fix -pd pred_data_adj_csem_mpi.txt -gd adj_csem_data_grad_mpi.txt -es esolution_adj_csem_mpi.bin -c ../docs/control_file_template
# #
# # SERIAL VALGRIND
# #
# valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=csem_valgrind_fwd_serial.txt ./TestSerial -f -m ../inputs/Naser_CSEM/rFile_Model -d ../inputs/Naser_CSEM/rFile_Data_fix -pd pred_data_fwd_csem_serial.txt -es esolution_fwd_csem_serial.bin -c ../docs/control_file_template
# #
# valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=csem_valgrind_adjt_serial.txt ./TestSerial -jt -m ../inputs/Naser_CSEM/rFile_Model -d ../inputs/Naser_CSEM/rFile_Data_fix -pd pred_data_adjt_csem_serial.txt -dm adjt_csem_sigma_serial.txt -es esolution_adjt_csem_serial.bin -c ../docs/control_file_template
# #
# valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=csem_valgrind_adj_serial.txt ./TestSerial -j -m ../inputs/Naser_CSEM/rFile_Model -pm adjt_csem_sigma_serial.txt -d ../inputs/Naser_CSEM/rFile_Data_fix -pd pred_data_adj_csem_serial.txt -gd adj_csem_data_grad_serial.txt -es esolution_adj_csem_serial.bin -c ../docs/control_file_template
# #
# MPI
#
mpirun -np 46 ./TestMPI -f -m ../inputs/Naser_CSEM/rFile_Model -d ../inputs/Naser_CSEM/rFile_Data_fix -pd pred_data_fwd_csem_mpi.txt -es esolution_fwd_csem_mpi.bin -c ../docs/control_file_template
#
mpirun -np 46 ./TestMPI -jt -m ../inputs/Naser_CSEM/rFile_Model -d ../inputs/Naser_CSEM/rFile_Data_fix -pd pred_data_adjt_csem_mpi.txt -dm adjt_csem_sigma_mpi.txt -es esolution_adjt_csem_mpi.bin -c ../docs/control_file_template
#
mpirun -np 46 ./TestMPI -j -m ../inputs/Naser_CSEM/rFile_Model -pm adjt_csem_sigma_mpi.txt -d ../inputs/Naser_CSEM/rFile_Data_fix -pd pred_data_adj_csem_mpi.txt -gd adj_csem_data_grad_mpi.txt -es esolution_adj_csem_mpi.bin -c ../docs/control_file_template
#
# SERIAL
#
#./TestSerial -f -m ../inputs/Naser_CSEM/rFile_Model -d ../inputs/Naser_CSEM/rFile_Data_fix -pd pred_data_fwd_csem_serial.txt -es esolution_fwd_csem_serial.bin -c ../docs/control_file_template
#
#./TestSerial -jt -m ../inputs/Naser_CSEM/rFile_Model -d ../inputs/Naser_CSEM/rFile_Data_fix -pd pred_data_adjt_csem_serial.txt -dm adjt_csem_sigma_serial.txt -es esolution_adjt_csem_serial.bin -c ../docs/control_file_template
#
#./TestSerial -j -m ../inputs/Naser_CSEM/rFile_Model -pm adjt_csem_sigma_serial.txt -d ../inputs/Naser_CSEM/rFile_Data_fix -pd pred_data_adj_csem_serial.txt -gd adj_csem_data_grad_serial.txt -es esolution_adj_csem_serial.bin -c ../docs/control_file_template
#
#mpirun -np 2 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=csem1_valgrind_fwd_mpi.txt ./TestMPI -f -m ../inputs/1st_Example/rFile_Model -d ../inputs/1st_Example/rFile_Data_CSEM -pd pred_data_fwd_csem1_mpi.txt -es esolution_fwd_csem1_mpi.bin -c ../docs/control_file_template
#
#mpirun -np 2 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=csem1_valgrind_adjt_mpi.txt ./TestMPI -jt -m ../inputs/1st_Example/rFile_Model -d ../inputs/1st_Example/rFile_Data_CSEM -pd pred_data_adjt_csem1_mpi.txt -dm adjt_csem1_sigma_mpi.txt -es esolution_adjt_csem1_mpi.bin -c ../docs/control_file_template
#
#mpirun -np 2 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=csem1_valgrind_adj_mpi.txt ./TestMPI -j -m ../inputs/1st_Example/rFile_Model -pm adjt_csem1_sigma_mpi.txt -d ../inputs/1st_Example/rFile_Data_CSEM -pd pred_data_adj_csem1_mpi.txt -es esolution_adj_csem1_mpi.bin -c ../docs/control_file_template