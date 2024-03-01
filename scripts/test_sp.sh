#!/bin/bash
#
########
# MT
########
#
# FWD
#
mpirun -np 13 --hostfile ../../hostfile_13.txt ./ModEM_MPI -f -m ../../modem_inputs/MT_5BLOCKS/m_VTI.rho -d ../../modem_inputs/MT_5BLOCKS/d0.dat -cf ../../modem_inputs/Benchmarking/oo_fwd_control_mf_sg.txt -pd MT_5BLOCKS_PRED_DATA_MF.dat
#
mpirun -np 13 --hostfile ../../hostfile_13.txt ./ModEM_MPI -f -m ../../modem_inputs/MT_5BLOCKS/m_VTI.rho -d ../../modem_inputs/MT_5BLOCKS/d0.dat -cf ../../modem_inputs/Benchmarking/oo_fwd_control_sp_sg.txt -pd MT_5BLOCKS_PRED_DATA_SP.dat
#
# JMULT
#
mpirun -np 13 --hostfile ../../hostfile_13.txt ./ModEM_MPI -j -m ../../modem_inputs/MT_5BLOCKS/m0_VTI.rho -pm ../../modem_inputs/MT_5BLOCKS/m_VTI.rho -d ../../modem_inputs/MT_5BLOCKS/d0.dat -cf ../../modem_inputs/Benchmarking/oo_fwd_control_mf_sg.txt -jm MT_5BLOCKS_JMHAT_MF.dat
#
mpirun -np 13 --hostfile ../../hostfile_13.txt ./ModEM_MPI -j -m ../../modem_inputs/MT_5BLOCKS/m0_VTI.rho -pm ../../modem_inputs/MT_5BLOCKS/m_VTI.rho -d ../../modem_inputs/MT_5BLOCKS/d0.dat -cf ../../modem_inputs/Benchmarking/oo_fwd_control_sp_sg.txt -jm MT_5BLOCKS_JMHAT_SP.dat
#
# JMULT_T
#
mpirun -np 13 --hostfile ../../hostfile_13.txt ./ModEM_MPI -jt -m ../../modem_inputs/MT_5BLOCKS/m0_VTI.rho -d ../../modem_inputs/MT_5BLOCKS/d_05.dat -cf ../../modem_inputs/Benchmarking/oo_fwd_control_mf_sg.txt -dm MT_5BLOCKS_DSIGMA_MF.rho
#
mpirun -np 13 --hostfile ../../hostfile_13.txt ./ModEM_MPI -jt -m ../../modem_inputs/MT_5BLOCKS/m0_VTI.rho -d ../../modem_inputs/MT_5BLOCKS/d_05.dat -cf ../../modem_inputs/Benchmarking/oo_fwd_control_sp_sg.txt -dm MT_5BLOCKS_DSIGMA_SP.rho
#
# INV
#
mpirun -np 13 --hostfile ../../hostfile_13.txt ./ModEM_MPI -i -m ../../modem_inputs/MT_5BLOCKS/m0_VTI.rho -d ../../modem_inputs/MT_5BLOCKS/d_05.dat -ci ../../modem_inputs/Benchmarking/oo_inv_control_nlcg -cf ../../modem_inputs/Benchmarking/oo_fwd_control_mf_sg.txt -o MT_5BLOCKS_INV_MF.rho
#
mpirun -np 13 --hostfile ../../hostfile_13.txt ./ModEM_MPI -i -m ../../modem_inputs/MT_5BLOCKS/m0_VTI.rho -d ../../modem_inputs/MT_5BLOCKS/d_05.dat -ci ../../modem_inputs/Benchmarking/oo_inv_control_nlcg -cf ../../modem_inputs/Benchmarking/oo_fwd_control_sp_sg.txt -o MT_5BLOCKS_INV_SP.rho
#
########
# CSEM
########
#
# FWD
#
mpirun -np 55 --hostfile ../../hostfile_55.txt ./ModEM_MPI -f -m ../../modem_inputs/CSEM_Small_model/Starting_model_Small_v_4_h.rho -d ../../modem_inputs/CSEM_Small_model/Data_on_Small_Grid_10%Err.dat -cf ../../modem_inputs/Benchmarking/oo_fwd_control_mf_sg.txt -pd CSEM_SMALL_PRED_DATA_MF.dat
#
mpirun -np 55 --hostfile ../../hostfile_55.txt ./ModEM_MPI -f -m ../../modem_inputs/CSEM_Small_model/Starting_model_Small_v_4_h.rho -d ../../modem_inputs/CSEM_Small_model/Data_on_Small_Grid_10%Err.dat -cf ../../modem_inputs/Benchmarking/oo_fwd_control_sp_sg.txt -pd CSEM_SMALL_PRED_DATA_SP.dat
#
# JMULT
#
mpirun -np 55 --hostfile ../../hostfile_55.txt ./ModEM_MPI -j -m ../../modem_inputs/CSEM_Small_model/Starting_model_Small_v_4_h.rho -pm ../../modem_inputs/CSEM_Small_model/Starting_model_Small_v_4_h.rho -d ../../modem_inputs/CSEM_Small_model/Data_on_Small_Grid_10%Err.dat -cf ../../modem_inputs/Benchmarking/oo_fwd_control_mf_sg.txt -jm CSEM_SMALL_JMHAT_MF.dat
#
mpirun -np 55 --hostfile ../../hostfile_55.txt ./ModEM_MPI -j -m ../../modem_inputs/CSEM_Small_model/Starting_model_Small_v_4_h.rho -pm ../../modem_inputs/CSEM_Small_model/Starting_model_Small_v_4_h.rho -d ../../modem_inputs/CSEM_Small_model/Data_on_Small_Grid_10%Err.dat -cf ../../modem_inputs/Benchmarking/oo_fwd_control_sp_sg.txt -jm CSEM_SMALL_JMHAT_SP.dat
#
# JMULT_T
#
mpirun -np 55 --hostfile ../../hostfile_55.txt ./ModEM_MPI -jt -m ../../modem_inputs/CSEM_Small_model/Starting_model_Small_v_4_h.rho -d ../../modem_inputs/CSEM_Small_model/Data_on_Small_Grid_10%Err.dat -cf ../../modem_inputs/Benchmarking/oo_fwd_control_mf_sg.txt -dm CSEM_SMALL_DSIGMA_MF.rho
#
mpirun -np 55 --hostfile ../../hostfile_55.txt ./ModEM_MPI -jt -m ../../modem_inputs/CSEM_Small_model/Starting_model_Small_v_4_h.rho -d ../../modem_inputs/CSEM_Small_model/Data_on_Small_Grid_10%Err.dat -cf ../../modem_inputs/Benchmarking/oo_fwd_control_sp_sg.txt -dm CSEM_SMALL_DSIGMA_SP.rho
#
# INV
#
mpirun -np 55 --hostfile ../../hostfile_55.txt ./ModEM_MPI -i -m ../../modem_inputs/CSEM_Small_model/Starting_model_Small_v_4_h.rho -d ../../modem_inputs/CSEM_Small_model/Data_on_Small_Grid_10%Err.dat -ci ../../modem_inputs/Benchmarking/oo_inv_control_nlcg -cf ../../modem_inputs/Benchmarking/oo_fwd_control_mf_sg.txt -o CSEM_SMALL_INV_MF.rho
#
mpirun -np 55 --hostfile ../../hostfile_55.txt ./ModEM_MPI -i -m ../../modem_inputs/CSEM_Small_model/Starting_model_Small_v_4_h.rho -d ../../modem_inputs/CSEM_Small_model/Data_on_Small_Grid_10%Err.dat -ci ../../modem_inputs/Benchmarking/oo_inv_control_nlcg -cf ../../modem_inputs/Benchmarking/oo_fwd_control_sp_sg.txt -o CSEM_SMALL_INV_SP.rho
#