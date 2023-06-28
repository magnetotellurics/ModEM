#!/bin/bash
#
# OO FWD
#
# SERIAL
#
#./ModEM_SERIAL -f -m ../../results/model_ISO_h10.mod -d ../../results/Data_Z_1tx_1rx.dat -cf fwd_ctrl_em1d -pd pred_mt1_iso_oo_serial.dat
#
#./ModEM_SERIAL -f -m ../../results/VTI_H10_V10/model_VTI_h10_v10.mod -d ../../results/Data_Z_1tx_1rx.dat -cf fwd_ctrl_em1d -pd pred_mt1_vti_oo_serial.dat
#
#./ModEM_SERIAL -f -m ../../results/model_ISO_h10.mod -d ../../results/Data_Ey_1tx_1rx.dat -cf fwd_ctrl_dip1d -pd pred_dip1d_oo_serial.dat
#
#./ModEM_SERIAL -f -m ../../results/model_ISO_h10.mod -d ../../results/Data_Ey_1tx_1rx.dat -cf fwd_ctrl_em1d -pd pred_em1d_iso_oo_serial.dat
#
#./ModEM_SERIAL -f -m ../../results/VTI_H10_V10/model_VTI_h10_v10.mod -d ../../results/Data_Ey_1tx_1rx.dat -cf fwd_ctrl_em1d -pd pred_em1d_vti_oo_serial.dat
#
# MPI
#
mpirun -np 22 ./ModEM_MPI -f -m ../../results/model_ISO_h10.mod -d ../../results/Data_Z_1tx_1rx.dat -cf fwd_ctrl_em1d -pd pred_mt_iso_oo_mpi.dat
#
mpirun -np 22 ./ModEM_MPI -f -m ../../results/VTI_H10_V10/model_VTI_h10_v10.mod -d ../../results/Data_Z_1tx_1rx.dat -cf fwd_ctrl_em1d -pd pred_mt_vti_oo_mpi.dat
#
mpirun -np 22 ./ModEM_MPI -f -m ../../results/model_ISO_h10.mod -d ../../results/Data_Ey_1tx_1rx.dat -cf fwd_ctrl_dip1d -pd pred_dip1d_oo_mpi.dat
#
mpirun -np 22 ./ModEM_MPI -f -m ../../results/model_ISO_h10.mod -d ../../results/Data_Ey_1tx_1rx.dat -cf fwd_ctrl_em1d -pd pred_em1d_iso_oo_mpi.dat
#
mpirun -np 22 ./ModEM_MPI -f -m ../../results/VTI_H10_V10/model_VTI_h10_v10.mod -d ../../results/Data_Ey_1tx_1rx.dat -cf fwd_ctrl_em1d -pd pred_em1d_vti_oo_mpi.dat
#