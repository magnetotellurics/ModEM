#!/bin/bash
#
# OO
#
# FWD
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./Mod3DMT_STD -F ../../modem-oo/docs/vti_results/model_ISO_h10.mod ../../modem-oo/docs/vti_results/Data_Z.dat pred_mt_iso_on.dat esol_mt_iso_on.bin ../../modem-oo/docs/vti_results/on_fwd_control_em1d
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./Mod3DMT_STD -F ../../modem-oo/docs/vti_results/model_VTI_h10_v10.mod ../../modem-oo/docs/vti_results/Data_Z.dat pred_mt_vti_on.dat esol_mt_vti_on.bin ../../modem-oo/docs/vti_results/on_fwd_control_em1d
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./Mod3DMT_STD -F ../../modem-oo/docs/vti_results/model_ISO_h10.mod ../../modem-oo/docs/vti_results/Data_Ey.dat pred_dip1d_on.dat esol_dip1d_on.bin ../../modem-oo/docs/vti_results/on_fwd_control_dip1d
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./Mod3DMT_STD -F ../../modem-oo/docs/vti_results/model_ISO_h10.mod ../../modem-oo/docs/vti_results/Data_Ey.dat pred_em1d_iso_on.dat esol_em1d_iso_on.bin ../../modem-oo/docs/vti_results/on_fwd_control_em1d
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./Mod3DMT_STD -F ../../modem-oo/docs/vti_results/model_VTI_h10_v10.mod ../../modem-oo/docs/vti_results/Data_Ey.dat pred_em1d_vti_on.dat esol_em1d_vti_on.bin ../../modem-oo/docs/vti_results/on_fwd_control_em1d
#
# JMULT
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./Mod3DMT_STD -M ../../modem-oo/docs/vti_results/model_ISO_h0.mod ../../modem-oo/docs/vti_results/model_ISO_h10.mod ../../modem-oo/docs/vti_results/Data_Z.dat jmhat_mt_iso_on.dat ../../modem-oo/docs/vti_results/on_fwd_control_em1d
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./Mod3DMT_STD -M ../../modem-oo/docs/vti_results/model_VTI_h0_v0.mod ../../modem-oo/docs/vti_results/model_VTI_h10_v10.mod ../../modem-oo/docs/vti_results/Data_Z.dat jmhat_mt_vti_on.dat ../../modem-oo/docs/vti_results/on_fwd_control_em1d
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./Mod3DMT_STD -M ../../modem-oo/docs/vti_results/model_ISO_h0.mod ../../modem-oo/docs/vti_results/model_ISO_h10.mod ../../modem-oo/docs/vti_results/Data_Ey.dat jmhat_dip1d_on.dat ../../modem-oo/docs/vti_results/on_fwd_control_dip1d
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./Mod3DMT_STD -M ../../modem-oo/docs/vti_results/model_ISO_h0.mod ../../modem-oo/docs/vti_results/model_ISO_h10.mod ../../modem-oo/docs/vti_results/Data_Ey.dat jmhat_em1d_iso_on.dat ../../modem-oo/docs/vti_results/on_fwd_control_em1d
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./Mod3DMT_STD -M ../../modem-oo/docs/vti_results/model_VTI_h0_v0.mod ../../modem-oo/docs/vti_results/model_VTI_h10_v10.mod ../../modem-oo/docs/vti_results/Data_Ey.dat jmhat_em1d_vti_on.dat ../../modem-oo/docs/vti_results/on_fwd_control_em1d
#
# JMULT
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./Mod3DMT_STD -T ../../modem-oo/docs/vti_results/model_ISO_h0.mod pred_mt_iso_on.dat dsigma_mt_iso_on.rho ../../modem-oo/docs/vti_results/on_fwd_control_em1d
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./Mod3DMT_STD -T ../../modem-oo/docs/vti_results/model_VTI_h0_v0.mod pred_mt_vti_on.dat dsigma_mt_vti_on.rho ../../modem-oo/docs/vti_results/on_fwd_control_em1d
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./Mod3DMT_STD -T ../../modem-oo/docs/vti_results/model_ISO_h0.mod pred_dip1d_on.dat dsigma_dip1d_on.rho ../../modem-oo/docs/vti_results/on_fwd_control_dip1d
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./Mod3DMT_STD -T ../../modem-oo/docs/vti_results/model_ISO_h0.mod pred_em1d_iso_on.dat dsigma_em1d_iso_on.rho ../../modem-oo/docs/vti_results/on_fwd_control_em1d
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./Mod3DMT_STD -T ../../modem-oo/docs/vti_results/model_VTI_h0_v0.mod pred_em1d_vti_on.dat dsigma_em1d_vti_on.rho ../../modem-oo/docs/vti_results/on_fwd_control_em1d
#
# NLCG
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./Mod3DMT_STD -I NLCG ../../modem-oo/docs/vti_results/model_ISO_h0.mod pred_mt_iso_on.dat ../../modem-oo/docs/vti_results/on_inv_control_iso ../../modem-oo/docs/vti_results/on_fwd_control_em1d -o mt_iso_on_nlcg
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./Mod3DMT_STD -I NLCG ../../modem-oo/docs/vti_results/model_VTI_h0_v0.mod pred_mt_vti_on.dat ../../modem-oo/docs/vti_results/on_inv_control_vti ../../modem-oo/docs/vti_results/on_fwd_control_em1d -o mt_vti_on_nlcg
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./Mod3DMT_STD -I NLCG ../../modem-oo/docs/vti_results/model_ISO_h0.mod pred_dip1d_on.dat ../../modem-oo/docs/vti_results/on_inv_control_iso ../../modem-oo/docs/vti_results/on_fwd_control_dip1d -o dip1d_on_nlcg
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./Mod3DMT_STD -I NLCG ../../modem-oo/docs/vti_results/model_ISO_h0.mod pred_em1d_iso_on.dat ../../modem-oo/docs/vti_results/on_inv_control_iso ../../modem-oo/docs/vti_results/on_fwd_control_em1d -o em1d_iso_on_nlcg
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./Mod3DMT_STD -I NLCG ../../modem-oo/docs/vti_results/model_VTI_h0_v0.mod pred_em1d_vti_on.dat ../../modem-oo/docs/vti_results/on_inv_control_vti ../../modem-oo/docs/vti_results/on_fwd_control_em1d -o em1d_vti_on_nlcg
#
# DCG
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./Mod3DMT_STD -I DCG ../../modem-oo/docs/vti_results/model_ISO_h0.mod pred_mt_iso_on.dat ../../modem-oo/docs/vti_results/on_inv_control_iso ../../modem-oo/docs/vti_results/on_fwd_control_em1d -o mt_iso_on_nlcg
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./Mod3DMT_STD -I DCG ../../modem-oo/docs/vti_results/model_VTI_h0_v0.mod pred_mt_vti_on.dat ../../modem-oo/docs/vti_results/on_inv_control_vti ../../modem-oo/docs/vti_results/on_fwd_control_em1d -o mt_vti_on_nlcg
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./Mod3DMT_STD -I DCG ../../modem-oo/docs/vti_results/model_ISO_h0.mod pred_dip1d_on.dat ../../modem-oo/docs/vti_results/on_inv_control_iso ../../modem-oo/docs/vti_results/on_fwd_control_dip1d -o dip1d_on_nlcg
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./Mod3DMT_STD -I DCG ../../modem-oo/docs/vti_results/model_ISO_h0.mod pred_em1d_iso_on.dat ../../modem-oo/docs/vti_results/on_inv_control_iso ../../modem-oo/docs/vti_results/on_fwd_control_em1d -o em1d_iso_on_nlcg
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./Mod3DMT_STD -I DCG ../../modem-oo/docs/vti_results/model_VTI_h0_v0.mod pred_em1d_vti_on.dat ../../modem-oo/docs/vti_results/on_inv_control_vti ../../modem-oo/docs/vti_results/on_fwd_control_em1d -o em1d_vti_on_nlcg
#
# MPI VALGRIND
#
#mpirun -np 22 --hostfile ../../hostfile_22.txt valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=valgrind_pred_mt_iso_on.txt ./Mod3DMT_STD -F ../../modem-oo/docs/vti_results/model_ISO_h10.mod ../../modem-oo/docs/vti_results/Data_Z.dat ../../modem-oo/docs/vti_results/on_fwd_control_em1d pred_mt_iso_on.dat
#