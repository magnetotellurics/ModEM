#!/bin/bash
#
# INPUTS
#
procs=55
m0_iso="../../../modem_inputs/V_4H/Starting_model.rho_h"
pm_iso="../../../modem_inputs/V_4H/True_model.rho_h"
m0_vti="../../../modem_inputs/V_4H/Starting_model.rho"
pm_vti="../../../modem_inputs/V_4H/True_model.rho"
cov="../../../modem_inputs/V_4H/Starting_model.cov"
data="../../../modem_inputs/V_4H/Data_file.dat"
ctrl_fwd="oo_fwd_control_em1d"
ctrl_nlcg="oo_inv_control_nlcg"
hostfile="--hostfile ../../../hostfile_$procs.txt"
#hostfile=""
#
################################################################
# FWD
################################################################
#
#mpirun -np $procs $hostfile ../Mod3DMT_STD -F $pm_iso $data pred_mt_iso_oo.dat esol_mt_iso_oo.bin $ctrl_fwd
#
mpirun -np $procs $hostfile ../Mod3DMT_STD -F $pm_vti $data pred_mt_vti_oo.dat esol_mt_vti_oo.bin $ctrl_fwd
#
################################################################
# JMULT
################################################################
#
#mpirun -np $procs $hostfile ../Mod3DMT_STD -M $m0_iso $pm_iso $data jmhat_mt_iso_oo.dat $ctrl_fwd
#
mpirun -np $procs $hostfile ../Mod3DMT_STD -M $m0_vti $pm_vti $data jmhat_mt_vti_oo.dat $ctrl_fwd
#
################################################################
# JMULT_T
################################################################
#
#mpirun -np $procs $hostfile ../Mod3DMT_STD -T $m0_iso $datam dsigma_mt_iso_oo.rho $ctrl_fwd
#
mpirun -np $procs $hostfile ../Mod3DMT_STD -T $m0_vti $data dsigma_mt_vti_oo.rho $ctrl_fwd
#
################################################################
# NLCG
################################################################
#
#mpirun -np $procs $hostfile ../Mod3DMT_STD -I NLCG $m0_iso data $ctrl_nlcg $ctrl_fwd $cov
#
mpirun -np $procs $hostfile ../Mod3DMT_STD -I NLCG $m0_vti $data $ctrl_nlcg $ctrl_fwd $cov
#
