#!/bin/bash
#
# OO FWD
#
# SERIAL
#
#./ModEM_SERIAL -f -m ../../results/model_ISO_h10.mod -d ../../results/Data_Z.dat -cf fwd_ctrl_em1d -pd pred_mt1_iso_oo_serial.dat
#
#./ModEM_SERIAL -f -m ../../results/VTI_H10_V10/model_VTI_h10_v10.mod -d ../../results/Data_Z.dat -cf fwd_ctrl_em1d -pd pred_mt1_vti_oo_serial.dat
#
#./ModEM_SERIAL -f -m ../../results/model_ISO_h10.mod -d ../../results/Data_Ey.dat -cf fwd_ctrl_dip1d -pd pred_dip1d_oo_serial.dat
#
#./ModEM_SERIAL -f -m ../../results/model_ISO_h10.mod -d ../../results/Data_Ey.dat -cf fwd_ctrl_em1d -pd pred_em1d_iso_oo_serial.dat
#
#./ModEM_SERIAL -f -m ../../results/VTI_H10_V10/model_VTI_h10_v10.mod -d ../../results/Data_Ey.dat -cf fwd_ctrl_em1d -pd pred_em1d_vti_oo_serial.dat
#
# MPI
#
mpirun -np 22 ./ModEM_MPI -f -m ../../results/model_ISO_h10.mod -d ../../results/Data_Z.dat -cf fwd_ctrl_em1d -pd pred_mt_iso_oo_mpi.dat
#
mpirun -np 22 ./ModEM_MPI -f -m ../../results/VTI_H10_V10/model_VTI_h10_v10.mod -d ../../results/Data_Z.dat -cf fwd_ctrl_em1d -pd pred_mt_vti_oo_mpi.dat
#
mpirun -np 22 ./ModEM_MPI -f -m ../../results/model_ISO_h10.mod -d ../../results/Data_Ey.dat -cf fwd_ctrl_dip1d -pd pred_dip1d_oo_mpi.dat
#
mpirun -np 22 ./ModEM_MPI -f -m ../../results/model_ISO_h10.mod -d ../../results/Data_Ey.dat -cf fwd_ctrl_em1d -pd pred_em1d_iso_oo_mpi.dat
#
mpirun -np 22 ./ModEM_MPI -f -m ../../results/VTI_H10_V10/model_VTI_h10_v10.mod -d ../../results/Data_Ey.dat -cf fwd_ctrl_em1d -pd pred_em1d_vti_oo_mpi.dat
#
mpirun -np 22 ./ModEM_MPI -j -m ../../results/model_ISO_h0.mod -pm ../../results/model_ISO_h10.mod -d ../../results/Data_Z.dat -cf fwd_ctrl_em1d -jm jmhat_mt_iso_oo_mpi.dat
#
mpirun -np 22 ./ModEM_MPI -j -m ../../results/model_ISO_h0.mod -pm ../../results/VTI_H10_V10/model_VTI_h10_v10.mod -d ../../results/Data_Z.dat -cf fwd_ctrl_em1d -jm jmhat_mt_vti_oo_mpi.dat
#
mpirun -np 22 ./ModEM_MPI -jt -m ../../results/model_ISO_h10.mod -d ../../results/Data_Z.dat -cf fwd_ctrl_em1d -dm sigma_mt_iso_oo_mpi.rho
#
mpirun -np 22 ./ModEM_MPI -jt -m ../../results/VTI_H10_V10/model_VTI_h10_v10.mod -d ../../results/Data_Z.dat -cf fwd_ctrl_em1d -dm sigma_mt_vti_oo_mpi.rho
#
mpirun -np 22 ./ModEM_MPI -i -m ../../results/model_ISO_h10.mod -d ../../results/Data_Z.dat -cf fwd_ctrl_em1d -ci inv_contrl_nlcg.txt -o nlcg_mt_iso_oo_mpi
#
mpirun -np 22 ./ModEM_MPI -i -m ../../results/VTI_H10_V10/model_VTI_h10_v10.mod -d ../../results/Data_Z.dat -cf fwd_ctrl_em1d -ci inv_contrl_nlcg.txt -o nlcg_mt_vti_oo_mpi
#
mpirun -np 22 ./ModEM_MPI -i -m ../../results/model_ISO_h10.mod -d ../../results/Data_Z.dat -cf fwd_ctrl_em1d -ci inv_contrl_dcg.txt -o nlcg_mt_iso_oo_mpi
#
mpirun -np 22 ./ModEM_MPI -i -m ../../results/VTI_H10_V10/model_VTI_h10_v10.mod -d ../../results/Data_Z.dat -cf fwd_ctrl_em1d -ci inv_contrl_dcg.txt -o nlcg_mt_vti_oo_mpi
#
# MPI VALGRIND
#
#mpirun -np 22 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=valgrind_pred_mt_iso_oo_mpi.txt ./ModEM_MPI -f -m ../../results/model_ISO_h10.mod -d ../../results/Data_Z.dat -cf fwd_ctrl_em1d -pd pred_mt_iso_oo_mpi.dat
#
#mpirun -np 22 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=valgrind_pred_mt_vti_oo_mpi.txt ./ModEM_MPI -f -m ../../results/VTI_H10_V10/model_VTI_h10_v10.mod -d ../../results/Data_Z.dat -cf fwd_ctrl_em1d -pd pred_mt_vti_oo_mpi.dat
#
#mpirun -np 22 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=valgrind_pred_dip1d_oo_mpi.txt ./ModEM_MPI -f -m ../../results/model_ISO_h10.mod -d ../../results/Data_Ey.dat -cf fwd_ctrl_dip1d -pd pred_dip1d_oo_mpi.dat
#
#mpirun -np 22 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=valgrind_pred_em1d_iso_oo_mpi.txt ./ModEM_MPI -f -m ../../results/model_ISO_h10.mod -d ../../results/Data_Ey.dat -cf fwd_ctrl_em1d -pd pred_em1d_iso_oo_mpi.dat
#
#mpirun -np 22 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=valgrind_pred_em1d_vti_oo_mpi.txt ./ModEM_MPI -f -m ../../results/VTI_H10_V10/model_VTI_h10_v10.mod -d ../../results/Data_Ey.dat -cf fwd_ctrl_em1d -pd pred_em1d_vti_oo_mpi.dat
#