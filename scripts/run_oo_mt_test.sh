#!/bin/bash
#
# MPI
#
./TestSerial -f -m ../inputs/1st_Example/rFile_Model_trim -d ../inputs/1st_Example/rFile_Data_MT -pd pred_data_fwd_mt16_mpi.txt -es esolution_fwd_mt16_mpi.bin -c ../docs/control_file_template
#
./TestSerial -jt -m ../inputs/1st_Example/rFile_Model_trim -d ../inputs/1st_Example/rFile_Data_MT -pd pred_data_adjt_mt16_mpi.txt -dm adjt_mt16_sigma_mpi.txt -es esolution_adjt_mt16_mpi.bin -c ../docs/control_file_template
#
./TestSerial -j -m ../inputs/1st_Example/rFile_Model_homo10ohm -pm ../inputs/1st_Example/rFile_Model_trim -d ../inputs/1st_Example/rFile_Data_MT -pd pred_data_adj_mt16_mpi.txt -gd adj_mt16_data_grad.txt -es esolution_adj_mt16_mpi.bin -c ../docs/control_file_template
#
# SERIAL
#
#valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=mt16_valgrind_fwd_serial.txt ./TestSerial -f -m ../inputs/1st_Example/rFile_Model_trim -d ../inputs/1st_Example/rFile_Data_MT -pd pred_data_fwd_mt16_serial.txt -es esolution_fwd_mt16_serial.bin -c ../docs/control_file_template
#
#valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=mt16_valgrind_adjt_serial.txt ./TestSerial -jt -m ../inputs/1st_Example/rFile_Model_trim -d ../inputs/1st_Example/rFile_Data_MT -pd pred_data_adjt_mt16_serial.txt -dm adjt_mt16_sigma_serial.txt -es esolution_adjt_mt16_serial.bin -c ../docs/control_file_template
#
#valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=mt16_valgrind_adj_serial.txt ./TestSerial -j -m ../inputs/1st_Example/rFile_Model_trim -pm adjt_mt16_sigma_serial.txt -d ../inputs/1st_Example/rFile_Data_MT -pd pred_data_adj_mt16_serial.txt -gd adj_mt16_data_grad.txt -es esolution_adj_mt16_serial.bin -c ../docs/control_file_template
#