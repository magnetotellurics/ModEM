#!/bin/bash
#
# OO
#
# FWD
#
mpirun -np 2 ./ModEM_MPI -f -m ../docs/vti_results/model_ISO_h10.mod -d ../docs/vti_results/Data_Z_1tx_1rx.dat -pd pred_mt_iso_oo.dat -es esol_mt_iso_oo.bin -cf ../docs/vti_results/oo_fwd_control_em1d
#
mpirun -np 2 ./ModEM_MPI -f -m ../docs/vti_results/VTI_H10_V10/model_VTI_h10_v10.mod -d ../docs/vti_results/Data_Z_1tx_1rx.dat -pd pred_mt_vti_oo.dat -es esol_mt_vti_oo.bin -cf ../docs/vti_results/oo_fwd_control_em1d
#
mpirun -np 2 ./ModEM_MPI -f -m ../docs/vti_results/model_ISO_h10.mod -d ../docs/vti_results/Data_Ey_1tx_1rx.dat -pd pred_dip1d_oo.dat -es esol_dip1d_oo.bin -cf ../docs/vti_results/oo_fwd_control_dip1d
#
mpirun -np 2 ./ModEM_MPI -f -m ../docs/vti_results/model_ISO_h10.mod -d ../docs/vti_results/Data_Ey_1tx_1rx.dat -pd pred_em1d_iso_oo.dat -es esol_em1d_iso_oo.bin -cf ../docs/vti_results/oo_fwd_control_em1d
#
mpirun -np 2 ./ModEM_MPI -f -m ../docs/vti_results/VTI_H10_V10/model_VTI_h10_v10.mod -d ../docs/vti_results/Data_Ey_1tx_1rx.dat -pd pred_em1d_vti_oo.dat -es esol_em1d_vti_oo.bin -cf ../docs/vti_results/oo_fwd_control_em1d
#
# JMULT
#
mpirun -np 2 ./ModEM_MPI -j -m ../docs/vti_results/model_ISO_h0.mod -pm ../docs/vti_results/model_ISO_h10.mod -d ../docs/vti_results/Data_Z_1tx_1rx.dat -jm jmhat_mt_iso_oo.dat -cf ../docs/vti_results/oo_fwd_control_em1d
#
mpirun -np 2 ./ModEM_MPI -j -m ../docs/vti_results/model_VTI_h0_v0.mod -pm ../docs/vti_results/VTI_H10_V10/model_VTI_h10_v10.mod -d ../docs/vti_results/Data_Z_1tx_1rx.dat -jm jmhat_mt_vti_oo.dat -cf ../docs/vti_results/oo_fwd_control_em1d
#
mpirun -np 2 ./ModEM_MPI -j -m ../docs/vti_results/model_ISO_h0.mod -pm ../docs/vti_results/model_ISO_h10.mod -d ../docs/vti_results/Data_Ey_1tx_1rx.dat -jm jmhat_dip1d_oo.dat -cf ../docs/vti_results/oo_fwd_control_dip1d
#
mpirun -np 2 ./ModEM_MPI -j -m ../docs/vti_results/model_ISO_h0.mod -pm ../docs/vti_results/model_ISO_h10.mod -d ../docs/vti_results/Data_Ey_1tx_1rx.dat -jm jmhat_em1d_iso_oo.dat -cf ../docs/vti_results/oo_fwd_control_em1d
#
mpirun -np 2 ./ModEM_MPI -j -m ../docs/vti_results/model_VTI_h0_v0.mod -pm ../docs/vti_results/VTI_H10_V10/model_VTI_h10_v10.mod -d ../docs/vti_results/Data_Ey_1tx_1rx.dat -jm jmhat_em1d_vti_oo.dat -cf ../docs/vti_results/oo_fwd_control_em1d
#
# JMULT
#
mpirun -np 2 ./ModEM_MPI -jt -m ../docs/vti_results/model_ISO_h0.mod -d pred_mt_iso_oo.dat -dm dsigma_mt_iso_oo.rho -cf ../docs/vti_results/oo_fwd_control_em1d
#
mpirun -np 2 ./ModEM_MPI -jt -m ../docs/vti_results/model_VTI_h0_v0.mod -d pred_mt_vti_oo.dat -dm dsigma_mt_vti_oo.rho -cf ../docs/vti_results/oo_fwd_control_em1d
#
mpirun -np 2 ./ModEM_MPI -jt -m ../docs/vti_results/model_ISO_h0.mod -d pred_dip1d_oo.dat -dm dsigma_dip1d_oo.rho -cf ../docs/vti_results/oo_fwd_control_dip1d
#
mpirun -np 2 ./ModEM_MPI -jt -m ../docs/vti_results/model_ISO_h0.mod -d pred_em1d_iso_oo.dat -dm dsigma_em1d_iso_oo.rho -cf ../docs/vti_results/oo_fwd_control_em1d
#
mpirun -np 2 ./ModEM_MPI -jt -m ../docs/vti_results/model_VTI_h0_v0.mod -d pred_em1d_vti_oo.dat -dm dsigma_em1d_vti_oo.rho -cf ../docs/vti_results/oo_fwd_control_em1d
#
# NLCG
#
mpirun -np 2 ./ModEM_MPI -i -m ../docs/vti_results/model_ISO_h0.mod -d pred_mt_iso_oo.dat -ci ../docs/vti_results/oo_inv_control_nlcg -cf ../docs/vti_results/oo_fwd_control_em1d -o mt_iso_oo_nlcg
#
mpirun -np 2 ./ModEM_MPI -i -m ../docs/vti_results/model_VTI_h0_v0.mod -d pred_mt_vti_oo.dat -ci ../docs/vti_results/oo_inv_control_nlcg -cf ../docs/vti_results/oo_fwd_control_em1d -o mt_vti_oo_nlcg
#
mpirun -np 2 ./ModEM_MPI -i -m ../docs/vti_results/model_ISO_h0.mod -d pred_dip1d_oo.dat -ci ../docs/vti_results/oo_inv_control_nlcg -cf ../docs/vti_results/oo_fwd_control_dip1d -o dip1d_oo_nlcg
#
mpirun -np 2 ./ModEM_MPI -i -m ../docs/vti_results/model_ISO_h0.mod -d pred_em1d_iso_oo.dat -ci ../docs/vti_results/oo_inv_control_nlcg -cf ../docs/vti_results/oo_fwd_control_em1d -o em1d_iso_oo_nlcg
#
mpirun -np 2 ./ModEM_MPI -i -m ../docs/vti_results/model_VTI_h0_v0.mod -d pred_em1d_vti_oo.dat -ci ../docs/vti_results/oo_inv_control_nlcg -cf ../docs/vti_results/oo_fwd_control_em1d -o em1d_vti_oo_nlcg
#
# DCG
#
mpirun -np 2 ./ModEM_MPI -i -m ../docs/vti_results/model_ISO_h0.mod -d pred_mt_iso_oo.dat -ci ../docs/vti_results/oo_inv_control_dcg -cf ../docs/vti_results/oo_fwd_control_em1d -o mt_iso_oo_nlcg
#
mpirun -np 2 ./ModEM_MPI -i -m ../docs/vti_results/model_VTI_h0_v0.mod -d pred_mt_vti_oo.dat -ci ../docs/vti_results/oo_inv_control_dcg -cf ../docs/vti_results/oo_fwd_control_em1d -o mt_vti_oo_nlcg
#
mpirun -np 2 ./ModEM_MPI -i -m ../docs/vti_results/model_ISO_h0.mod -d pred_dip1d_oo.dat -ci ../docs/vti_results/oo_inv_control_dcg -cf ../docs/vti_results/oo_fwd_control_dip1d -o dip1d_oo_nlcg
#
mpirun -np 2 ./ModEM_MPI -i -m ../docs/vti_results/model_ISO_h0.mod -d pred_em1d_iso_oo.dat -ci ../docs/vti_results/oo_inv_control_dcg -cf ../docs/vti_results/oo_fwd_control_em1d -o em1d_iso_oo_nlcg
#
mpirun -np 2 ./ModEM_MPI -i -m ../docs/vti_results/model_VTI_h0_v0.mod -d pred_em1d_vti_oo.dat -ci ../docs/vti_results/oo_inv_control_dcg -cf ../docs/vti_results/oo_fwd_control_em1d -o em1d_vti_oo_nlcg
#
# MPI VALGRIND
#
#mpirun -np 2 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=valgrind_pred_mt_iso_oo.txt ./ModEM_MPI -f -m ../docs/vti_results/model_ISO_h10.mod -d ../docs/vti_results/Data_Z_1tx_1rx.dat -cf ../docs/vti_results/oo_fwd_control_em1d -pd pred_mt_iso_oo.dat
#