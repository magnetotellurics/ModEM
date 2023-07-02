#!/bin/bash
#
################################################################
# FWD
################################################################
#
mpirun -np 2 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=pred_mt_iso_on.vg ./Mod3DMT_STD -F ../../modem-oo/docs/vti_results/model_ISO_h10.mod ../../modem-oo/docs/vti_results/Data_Z_1tx_1rx.dat pred_mt_iso_on.dat esol_mt_iso_on.bin on_fwd_control_em1d
#
mpirun -np 2 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=pred_mt_vti_on.vg ./Mod3DMT_STD -F ../../modem-oo/docs/vti_results/model_VTI_h10_v10.mod ../../modem-oo/docs/vti_results/Data_Z_1tx_1rx.dat pred_mt_vti_on.dat esol_mt_vti_on.bin on_fwd_control_em1d
#
mpirun -np 2 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=pred_dip1d_on.vg ./Mod3DMT_STD -F ../../modem-oo/docs/vti_results/model_ISO_h10.mod ../../modem-oo/docs/vti_results/Data_Ey_1tx_1rx.dat pred_dip1d_on.dat esol_dip1d_on.bin on_fwd_control_dip1d
#
mpirun -np 2 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=pred_em1d_iso_on.vg ./Mod3DMT_STD -F ../../modem-oo/docs/vti_results/model_ISO_h10.mod ../../modem-oo/docs/vti_results/Data_Ey_1tx_1rx.dat pred_em1d_iso_on.dat esol_em1d_iso_on.bin on_fwd_control_em1d
#
mpirun -np 2 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=pred_em1d_vti_on.vg ./Mod3DMT_STD -F ../../modem-oo/docs/vti_results/model_VTI_h10_v10.mod ../../modem-oo/docs/vti_results/Data_Ey_1tx_1rx.dat pred_em1d_vti_on.dat esol_em1d_vti_on.bin on_fwd_control_em1d
#
################################################################
# JMULT
################################################################
#
mpirun -np 2 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=jmhat_mt_iso_on.vg ./Mod3DMT_STD -M ../../modem-oo/docs/vti_results/model_ISO_h0.mod -pm ../../modem-oo/docs/vti_results/model_ISO_h10.mod ../../modem-oo/docs/vti_results/Data_Z_1tx_1rx.dat jmhat_mt_iso_on.dat on_fwd_control_em1d
#
mpirun -np 2 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=jmhat_mt_vti_on.vg ./Mod3DMT_STD -M ../../modem-oo/docs/vti_results/model_VTI_h0_v0.mod -pm ../../modem-oo/docs/vti_results/model_VTI_h10_v10.mod ../../modem-oo/docs/vti_results/Data_Z_1tx_1rx.dat jmhat_mt_vti_on.dat on_fwd_control_em1d
#
mpirun -np 2 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=jmhat_dip1d_on.vg ./Mod3DMT_STD -M ../../modem-oo/docs/vti_results/model_ISO_h0.mod -pm ../../modem-oo/docs/vti_results/model_ISO_h10.mod ../../modem-oo/docs/vti_results/Data_Ey_1tx_1rx.dat jmhat_dip1d_on.dat on_fwd_control_dip1d
#
mpirun -np 2 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=jmhat_em1d_iso_on.vg ./Mod3DMT_STD -M ../../modem-oo/docs/vti_results/model_ISO_h0.mod -pm ../../modem-oo/docs/vti_results/model_ISO_h10.mod ../../modem-oo/docs/vti_results/Data_Ey_1tx_1rx.dat jmhat_em1d_iso_on.dat on_fwd_control_em1d
#
mpirun -np 2 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=jmhat_em1d_vti_on.vg ./Mod3DMT_STD -M ../../modem-oo/docs/vti_results/model_VTI_h0_v0.mod -pm ../../modem-oo/docs/vti_results/model_VTI_h10_v10.mod ../../modem-oo/docs/vti_results/Data_Ey_1tx_1rx.dat jmhat_em1d_vti_on.dat on_fwd_control_em1d
#
################################################################
# JMULT_T
################################################################
#
mpirun -np 2 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=dsigma_mt_iso_on.vg ./Mod3DMT_STD -T ../../modem-oo/docs/vti_results/model_ISO_h0.mod pred_mt_iso_on.dat dsigma_mt_iso_on.rho on_fwd_control_em1d
#
mpirun -np 2 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=dsigma_mt_vti_on.vg ./Mod3DMT_STD -T ../../modem-oo/docs/vti_results/model_VTI_h0_v0.mod pred_mt_vti_on.dat dsigma_mt_vti_on.rho on_fwd_control_em1d
#
mpirun -np 2 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=dsigma_dip1d_on.vg ./Mod3DMT_STD -T ../../modem-oo/docs/vti_results/model_ISO_h0.mod pred_dip1d_on.dat dsigma_dip1d_on.rho on_fwd_control_dip1d
#
mpirun -np 2 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=dsigma_em1d_iso_on.vg ./Mod3DMT_STD -T ../../modem-oo/docs/vti_results/model_ISO_h0.mod pred_em1d_iso_on.dat dsigma_em1d_iso_on.rho on_fwd_control_em1d
#
mpirun -np 2 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=dsigma_em1d_vti_on.vg ./Mod3DMT_STD -T ../../modem-oo/docs/vti_results/model_VTI_h0_v0.mod pred_em1d_vti_on.dat dsigma_em1d_vti_on.rho on_fwd_control_em1d
#
################################################################
# NLCG
################################################################
#
mpirun -np 2 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=mt_iso_on_nlcg.vg ./Mod3DMT_STD -I NLCG ../../modem-oo/docs/vti_results/model_ISO_h0.mod pred_mt_iso_on.dat on_inv_control on_fwd_control_em1d
#
mpirun -np 2 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=mt_vti_on_nlcg.vg ./Mod3DMT_STD -I NLCG ../../modem-oo/docs/vti_results/model_VTI_h0_v0.mod pred_mt_vti_on.dat on_inv_control on_fwd_control_em1d
#
mpirun -np 2 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=dip1d_on_nlcg.vg ./Mod3DMT_STD -I NLCG ../../modem-oo/docs/vti_results/model_ISO_h0.mod pred_dip1d_on.dat on_inv_control on_fwd_control_dip1d
#
mpirun -np 2 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=em1d_iso_on_nlcg.vg ./Mod3DMT_STD -I NLCG ../../modem-oo/docs/vti_results/model_ISO_h0.mod pred_em1d_iso_on.dat on_inv_control on_fwd_control_em1d
#
mpirun -np 2 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=em1d_vti_on_nlcg.vg ./Mod3DMT_STD -I NLCG ../../modem-oo/docs/vti_results/model_VTI_h0_v0.mod pred_em1d_vti_on.dat on_inv_control on_fwd_control_em1d
#
################################################################
# DCG
################################################################
#
mpirun -np 2 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=mt_iso_on_dcg.vg ./Mod3DMT_STD -I DCG ../../modem-oo/docs/vti_results/model_ISO_h0.mod pred_mt_iso_on.dat on_inv_control on_fwd_control_em1d
#
mpirun -np 2 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=mt_vti_on_dcg.vg ./Mod3DMT_STD -I DCG ../../modem-oo/docs/vti_results/model_VTI_h0_v0.mod pred_mt_vti_on.dat on_inv_control on_fwd_control_em1d
#
mpirun -np 2 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=dip1d_on_dcg.vg ./Mod3DMT_STD -I DCG ../../modem-oo/docs/vti_results/model_ISO_h0.mod pred_dip1d_on.dat on_inv_control on_fwd_control_dip1d
#
mpirun -np 2 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=em1d_iso_on_dcg.vg ./Mod3DMT_STD -I DCG ../../modem-oo/docs/vti_results/model_ISO_h0.mod pred_em1d_iso_on.dat on_inv_control on_fwd_control_em1d
#
mpirun -np 2 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=em1d_vti_on_dcg.vg ./Mod3DMT_STD -I DCG ../../modem-oo/docs/vti_results/model_VTI_h0_v0.mod pred_em1d_vti_on.dat on_inv_control on_fwd_control_em1d
#