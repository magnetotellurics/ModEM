#!/bin/bash
#
# FWD
# #
# mpirun -np 55 --hostfile ../../hostfile_55.txt ./ModEM_MPI -f -m ../../modem_inputs/CSEM_Small_model/Starting_model_Small_v_4_h.rho_h -d ../../modem_inputs/CSEM_Small_model/Data_on_Small_Grid_10%Err.dat -cf it_mf_control -dm pred_it_mf.dat
# #
# mpirun -np 55 --hostfile ../../hostfile_55.txt ./ModEM_MPI -f -m ../../modem_inputs/CSEM_Small_model/Starting_model_Small_v_4_h.rho_h -d ../../modem_inputs/CSEM_Small_model/Data_on_Small_Grid_10%Err.dat -cf it_sp_control -dm pred_it_sp.dat
# #
# mpirun -np 55 --hostfile ../../hostfile_55.txt ./ModEM_MPI -f -m ../../modem_inputs/CSEM_Small_model/Starting_model_Small_v_4_h.rho_h -d ../../modem_inputs/CSEM_Small_model/Data_on_Small_Grid_10%Err.dat -cf dc_mf_control -dm pred_dc_mf.dat
# #
# mpirun -np 55 --hostfile ../../hostfile_55.txt ./ModEM_MPI -f -m ../../modem_inputs/CSEM_Small_model/Starting_model_Small_v_4_h.rho_h -d ../../modem_inputs/CSEM_Small_model/Data_on_Small_Grid_10%Err.dat -cf dc_sp_control -dm pred_dc_sp.dat
# #
# # JMULT
# #
# mpirun -np 55 --hostfile ../../hostfile_55.txt ./ModEM_MPI -j -m ../../modem_inputs/CSEM_Small_model/Starting_model_Small_v_4_h.rho_h -pm ../../modem_inputs/CSEM_Small_model/Starting_model_Small_v_4_h.rho_h -d ../../modem_inputs/CSEM_Small_model/Data_on_Small_Grid_10%Err.dat -cf it_mf_control -jm jmhat_it_mf.dat
# #
# mpirun -np 55 --hostfile ../../hostfile_55.txt ./ModEM_MPI -j -m ../../modem_inputs/CSEM_Small_model/Starting_model_Small_v_4_h.rho_h -pm ../../modem_inputs/CSEM_Small_model/Starting_model_Small_v_4_h.rho_h -d ../../modem_inputs/CSEM_Small_model/Data_on_Small_Grid_10%Err.dat -cf it_sp_control -jm jmhat_it_sp.dat
# #
# mpirun -np 55 --hostfile ../../hostfile_55.txt ./ModEM_MPI -j -m ../../modem_inputs/CSEM_Small_model/Starting_model_Small_v_4_h.rho_h -pm ../../modem_inputs/CSEM_Small_model/Starting_model_Small_v_4_h.rho_h -d ../../modem_inputs/CSEM_Small_model/Data_on_Small_Grid_10%Err.dat -cf dc_mf_control -jm jmhat_dc_mf.dat
# #
# mpirun -np 55 --hostfile ../../hostfile_55.txt ./ModEM_MPI -j -m ../../modem_inputs/CSEM_Small_model/Starting_model_Small_v_4_h.rho_h -pm ../../modem_inputs/CSEM_Small_model/Starting_model_Small_v_4_h.rho_h -d ../../modem_inputs/CSEM_Small_model/Data_on_Small_Grid_10%Err.dat -cf dc_sp_control -jm jmhat_dc_sp.dat
# #
# # JMULT_T
# #
# mpirun -np 55 --hostfile ../../hostfile_55.txt ./ModEM_MPI -jt -m ../../modem_inputs/CSEM_Small_model/Starting_model_Small_v_4_h.rho_h -d ../../modem_inputs/CSEM_Small_model/Data_on_Small_Grid_10%Err.dat -cf it_mf_control -dm dsigma_it_mf.rho
# #
# mpirun -np 55 --hostfile ../../hostfile_55.txt ./ModEM_MPI -jt -m ../../modem_inputs/CSEM_Small_model/Starting_model_Small_v_4_h.rho_h -d ../../modem_inputs/CSEM_Small_model/Data_on_Small_Grid_10%Err.dat -cf it_sp_control -dm dsigma_it_sp.rho
# #
# mpirun -np 55 --hostfile ../../hostfile_55.txt ./ModEM_MPI -jt -m ../../modem_inputs/CSEM_Small_model/Starting_model_Small_v_4_h.rho_h -d ../../modem_inputs/CSEM_Small_model/Data_on_Small_Grid_10%Err.dat -cf dc_mf_control -dm dsigma_dc_mf.rho
# #
# mpirun -np 55 --hostfile ../../hostfile_55.txt ./ModEM_MPI -jt -m ../../modem_inputs/CSEM_Small_model/Starting_model_Small_v_4_h.rho_h -d ../../modem_inputs/CSEM_Small_model/Data_on_Small_Grid_10%Err.dat -cf dc_sp_control -dm dsigma_dc_sp.rho
# #
# INV
#
mpirun -np 55 --hostfile ../../hostfile_55.txt ./ModEM_MPI -i -m ../../modem_inputs/CSEM_Small_model/Starting_model_Small_v_4_h.rho_h -d ../../modem_inputs/CSEM_Small_model/Data_on_Small_Grid_10%Err.dat -cf it_mf_control -o it_mf
#
mpirun -np 55 --hostfile ../../hostfile_55.txt ./ModEM_MPI -i -m ../../modem_inputs/CSEM_Small_model/Starting_model_Small_v_4_h.rho_h -d ../../modem_inputs/CSEM_Small_model/Data_on_Small_Grid_10%Err.dat -cf it_sp_control -o it_sp
#
mpirun -np 55 --hostfile ../../hostfile_55.txt ./ModEM_MPI -i -m ../../modem_inputs/CSEM_Small_model/Starting_model_Small_v_4_h.rho_h -d ../../modem_inputs/CSEM_Small_model/Data_on_Small_Grid_10%Err.dat -cf dc_mf_control -o dc_mf
#
mpirun -np 55 --hostfile ../../hostfile_55.txt ./ModEM_MPI -i -m ../../modem_inputs/CSEM_Small_model/Starting_model_Small_v_4_h.rho_h -d ../../modem_inputs/CSEM_Small_model/Data_on_Small_Grid_10%Err.dat -cf dc_sp_control -o dc_sp
#