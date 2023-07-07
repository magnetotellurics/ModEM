#!/bin/bash
#
# INPUTS
#
procs=2
m0_iso="model_ISO_h0.mod"
pm_iso="model_ISO_h10.mod"
m0_vti="model_VTI_h0_v0.mod"
pm_vti="model_VTI_h10_v10.mod"
data_mt="Data_Z_1tx_1rx.dat"
data_csem="Data_Ey_Z_1tx_1rx.dat"
ctrl_dip1d="on_fwd_control_dip1d"
ctrl_em1d="on_fwd_control_em1d"
ctrl_inv="on_inv_control"
#hostfile="--hostfile ../../hostfile_$procs.txt"
hostfile=""
#
# FWD
#
mpirun -np $procs $hostfile valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=pred_mt_iso_on.vg ./Mod3DMT_STD -F $pm_iso $data_mt pred_mt_iso_on.dat esol_mt_iso_on.bin $ctrl_em1d
#
mpirun -np $procs $hostfile valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=pred_mt_vti_on.vg ./Mod3DMT_STD -F $pm_vti $data_mt pred_mt_vti_on.dat esol_mt_vti_on.bin $ctrl_em1d
#
mpirun -np $procs $hostfile valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=pred_dip1d_on.vg ./Mod3DMT_STD -F $pm_iso $data_csem pred_dip1d_on.dat esol_dip1d_on.bin $ctrl_dip1d
#
mpirun -np $procs $hostfile valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=pred_em1d_iso_on.vg ./Mod3DMT_STD -F $pm_iso $data_csem pred_em1d_iso_on.dat esol_em1d_iso_on.bin $ctrl_em1d
#
mpirun -np $procs $hostfile valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=pred_em1d_vti_on.vg ./Mod3DMT_STD -F $pm_vti $data_csem pred_em1d_vti_on.dat esol_em1d_vti_on.bin $ctrl_em1d
#
# JMULT
#
mpirun -np $procs $hostfile valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=jmhat_mt_iso_on.vg ./Mod3DMT_STD -M $m0_iso $pm_iso $data_mt jmhat_mt_iso_on.dat $ctrl_em1d
#
mpirun -np $procs $hostfile valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=jmhat_mt_vti_on.vg ./Mod3DMT_STD -M $m0_vti $pm_vti $data_mt jmhat_mt_vti_on.dat $ctrl_em1d
#
mpirun -np $procs $hostfile valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=jmhat_dip1d_on.vg ./Mod3DMT_STD -M $m0_iso $pm_iso $data_csem jmhat_dip1d_on.dat $ctrl_dip1d
#
mpirun -np $procs $hostfile valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=jmhat_em1d_iso_on.vg ./Mod3DMT_STD -M $m0_iso $pm_iso $data_csem jmhat_em1d_iso_on.dat $ctrl_em1d
#
mpirun -np $procs $hostfile valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=jmhat_em1d_vti_on.vg ./Mod3DMT_STD -M $m0_vti $pm_vti $data_csem jmhat_em1d_vti_on.dat $ctrl_em1d
#
# JMULT_T
#
mpirun -np $procs $hostfile valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=dsigma_mt_iso_on.vg ./Mod3DMT_STD -T $m0_iso pred_mt_iso_on.dat dsigma_mt_iso_on.rho $ctrl_em1d
#
mpirun -np $procs $hostfile valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=dsigma_mt_vti_on.vg ./Mod3DMT_STD -T $m0_vti pred_mt_vti_on.dat dsigma_mt_vti_on.rho $ctrl_em1d
#
mpirun -np $procs $hostfile valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=dsigma_dip1d_on.vg ./Mod3DMT_STD -T $m0_iso pred_dip1d_on.dat dsigma_dip1d_on.rho $ctrl_dip1d
#
mpirun -np $procs $hostfile valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=dsigma_em1d_iso_on.vg ./Mod3DMT_STD -T $m0_iso pred_em1d_iso_on.dat dsigma_em1d_iso_on.rho $ctrl_em1d
#
mpirun -np $procs $hostfile valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=pred_mt_iso_on.vg ./Mod3DMT_STD -T $m0_vti pred_em1d_vti_on.dat dsigma_em1d_vti_on.rho $ctrl_em1d
#
# NLCG
#
mpirun -np $procs $hostfile valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=mt_iso_on.vg ./Mod3DMT_STD -I NLCG $m0_iso pred_mt_iso_on.dat $ctrl_inv $ctrl_em1d
#
#mpirun -np $procs $hostfile valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=pred_mt_iso_on.vg ./Mod3DMT_STD -I NLCG $m0_vti pred_mt_vti_on.dat $ctrl_inv $ctrl_em1d
#
#mpirun -np $procs $hostfile valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=pred_mt_iso_on.vg ./Mod3DMT_STD -I NLCG $m0_iso pred_dip1d_on.dat $ctrl_inv $ctrl_dip1d
#
#mpirun -np $procs $hostfile valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=pred_mt_iso_on.vg ./Mod3DMT_STD -I NLCG $m0_iso pred_em1d_iso_on.dat $ctrl_inv $ctrl_em1d
#
#mpirun -np $procs $hostfile valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=pred_mt_iso_on.vg ./Mod3DMT_STD -I NLCG $m0_vti pred_em1d_vti_on.dat $ctrl_inv $ctrl_em1d
#
# DCG
#
#mpirun -np $procs $hostfile valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=pred_mt_iso_on.vg ./Mod3DMT_STD -I DCG $m0_iso pred_mt_iso_on.dat $ctrl_inv $ctrl_em1d -o mt_iso_on_nlcg
#
#mpirun -np $procs $hostfile valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=pred_mt_iso_on.vg ./Mod3DMT_STD -I DCG $m0_vti pred_mt_vti_on.dat $ctrl_inv $ctrl_em1d -o mt_vti_on_nlcg
#
#mpirun -np $procs $hostfile valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=pred_mt_iso_on.vg ./Mod3DMT_STD -I DCG $m0_iso pred_dip1d_on.dat $ctrl_inv $ctrl_dip1d -o dip1d_on_nlcg
#
#mpirun -np $procs $hostfile valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=pred_mt_iso_on.vg ./Mod3DMT_STD -I DCG $m0_iso pred_em1d_iso_on.dat $ctrl_inv $ctrl_em1d -o em1d_iso_on_nlcg
#
#mpirun -np $procs $hostfile valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=pred_mt_iso_on.vg ./Mod3DMT_STD -I DCG $m0_vti pred_em1d_vti_on.dat $ctrl_inv $ctrl_em1d -o em1d_vti_on_nlcg
#
# MPI VALGRIND
#
#mpirun -np $procs $hostfile valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=valgrind_pred_mt_iso_on.txt ./Mod3DMT_STD -F $pm_iso $data_mt $ctrl_em1d pred_mt_iso_on.dat
#