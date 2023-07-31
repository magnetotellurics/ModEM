#!/bin/bash
#
# FWD
#
mpirun -np 5 ./ModEm_MPI -f -m ../../modem_inputs/Benchmarking/rFile_Model -d ../../modem_inputs/Benchmarking/rFile_Data -cf it_mf_control -dm pred_it_mf.dat
#
mpirun -np 5 ./ModEm_MPI -f -m ../../modem_inputs/Benchmarking/rFile_Model -d ../../modem_inputs/Benchmarking/rFile_Data -cf it_sp_control -dm pred_it_sp.dat
#
mpirun -np 5 ./ModEm_MPI -f -m ../../modem_inputs/Benchmarking/rFile_Model -d ../../modem_inputs/Benchmarking/rFile_Data -cf dc_mf_control -dm pred_dc_mf.dat
#
mpirun -np 5 ./ModEm_MPI -f -m ../../modem_inputs/Benchmarking/rFile_Model -d ../../modem_inputs/Benchmarking/rFile_Data -cf dc_sp_control -dm pred_dc_sp.dat
#
# JMULT
#
mpirun -np 5 ./ModEm_MPI -j -m ../../modem_inputs/Benchmarking/rFile_Model -pm ../../modem_inputs/Benchmarking/rFile_Model -d ../../modem_inputs/Benchmarking/rFile_Data -cf it_mf_control -jm jmhat_it_mf.dat
#
mpirun -np 5 ./ModEm_MPI -j -m ../../modem_inputs/Benchmarking/rFile_Model -pm ../../modem_inputs/Benchmarking/rFile_Model -d ../../modem_inputs/Benchmarking/rFile_Data -cf it_sp_control -jm jmhat_it_sp.dat
#
mpirun -np 5 ./ModEm_MPI -j -m ../../modem_inputs/Benchmarking/rFile_Model -pm ../../modem_inputs/Benchmarking/rFile_Model -d ../../modem_inputs/Benchmarking/rFile_Data -cf dc_mf_control -jm jmhat_dc_mf.dat
#
mpirun -np 5 ./ModEm_MPI -j -m ../../modem_inputs/Benchmarking/rFile_Model -pm ../../modem_inputs/Benchmarking/rFile_Model -d ../../modem_inputs/Benchmarking/rFile_Data -cf dc_sp_control -jm jmhat_dc_sp.dat
#
# JMULT_T
#
mpirun -np 5 ./ModEm_MPI -jt -m ../../modem_inputs/Benchmarking/rFile_Model -d ../../modem_inputs/Benchmarking/rFile_Data -cf it_mf_control -dm dsigma_it_mf.rho
#
mpirun -np 5 ./ModEm_MPI -jt -m ../../modem_inputs/Benchmarking/rFile_Model -d ../../modem_inputs/Benchmarking/rFile_Data -cf it_sp_control -dm dsigma_it_sp.rho
#
mpirun -np 5 ./ModEm_MPI -jt -m ../../modem_inputs/Benchmarking/rFile_Model -d ../../modem_inputs/Benchmarking/rFile_Data -cf dc_mf_control -dm dsigma_dc_mf.rho
#
mpirun -np 5 ./ModEm_MPI -jt -m ../../modem_inputs/Benchmarking/rFile_Model -d ../../modem_inputs/Benchmarking/rFile_Data -cf dc_sp_control -dm dsigma_dc_sp.rho
#
#INV
#
mpirun -np 5 ./ModEm_MPI -i -m ../../modem_inputs/Benchmarking/rFile_Model -d ../../modem_inputs/Benchmarking/rFile_Data -cf it_mf_control -o it_mf
#
mpirun -np 5 ./ModEm_MPI -i -m ../../modem_inputs/Benchmarking/rFile_Model -d ../../modem_inputs/Benchmarking/rFile_Data -cf it_sp_control -o it_sp
#
mpirun -np 5 ./ModEm_MPI -i -m ../../modem_inputs/Benchmarking/rFile_Model -d ../../modem_inputs/Benchmarking/rFile_Data -cf dc_mf_control -o dc_mf
#
mpirun -np 5 ./ModEm_MPI -i -m ../../modem_inputs/Benchmarking/rFile_Model -d ../../modem_inputs/Benchmarking/rFile_Data -cf dc_sp_control -o dc_sp
#