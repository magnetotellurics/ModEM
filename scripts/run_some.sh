#!/bin/bash
#
mpirun -np 5 ./TestMPI -f -m ../inputs/Benchmarking/rFile_Model -d ../inputs/Benchmarking/rFile_Data -pd mpi_mt4_fwd_pred.txt -es mpi_mt4_fwd_esolution.bin -c ../inputs/oo_lowtol_control_file
#
mpirun -np 5 ./TestMPI -jt -m ../inputs/Benchmarking/rFile_Model -d ../inputs/Benchmarking/rFile_Data -pd mpi_mt4_adjt_pred.txt -es mpi_adjt_esolution.bin -dm mpi_adjt_dsigma.txt -c ../inputs/oo_lowtol_control_file
#
mpirun -np 5 ./TestMPI -j -m ../inputs/Benchmarking/rFile_Model -pm input_dsigma.txt -d ../inputs/Benchmarking/rFile_Data -pd mpi_mt4_adj_pred.txt -es mpi_adj_esolution.bin -gd mpi_adj_grad_data.txt -c ../inputs/oo_lowtol_control_file
#
./TestSerial -f -m ../inputs/Benchmarking/rFile_Model -d ../inputs/Benchmarking/rFile_Data -pd serial_mt4_fwd_pred.txt -es serial_mt4_fwd_esolution.bin -c ../inputs/oo_lowtol_control_file
#
./TestSerial -jt -m ../inputs/Benchmarking/rFile_Model -d ../inputs/Benchmarking/rFile_Data -pd serial_mt4_adjt_pred.txt -es serial_adjt_esolution.bin -dm input_dsigma.txt -c ../inputs/oo_lowtol_control_file
#
./TestSerial -j -m ../inputs/Benchmarking/rFile_Model -pm input_dsigma.txt -d ../inputs/Benchmarking/rFile_Data -pd serial_mt4_adj_pred.txt -es serial_adj_esolution.bin -gd serial_adj_grad_data.txt -c ../inputs/oo_lowtol_control_file
#