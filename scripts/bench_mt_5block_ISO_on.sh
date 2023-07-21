#!/bin/bash
#
# INPUTS
#
procs=25
m0_iso="../../../modem_inputs/MT_5BLOCKS/m0.rho"
pm_iso="../../../modem_inputs/MT_5BLOCKS/m.rho"
m0_vti="../../../modem_inputs/MT_5BLOCKS/m0_VTI.rho"
pm_vti="../../../modem_inputs/MT_5BLOCKS/m_VTI.rho"
data_mt="../../../modem_inputs/MT_5BLOCKS/d_05.dat"
ctrl_fwd="on_fwd_control_em1d"
ctrl_nlcg="on_inv_control_nlcg"
ctrl_dcg="on_inv_control_dcg"
hostfile="--hostfile ../../../hostfile_$procs.txt"
#hostfile=""
#
################################################################
# FWD
################################################################
#
mpirun -np $procs $hostfile ../Mod3DMT_STD -F $pm_iso $data_mt pred_mt_iso_on.dat esol_mt_iso_on.bin $ctrl_fwd
#
#mpirun -np $procs $hostfile ../Mod3DMT_STD -F $pm_vti $data_mt pred_mt_vti_on.dat esol_mt_vti_on.bin $ctrl_fwd
#
################################################################
# JMULT
################################################################
#
mpirun -np $procs $hostfile ../Mod3DMT_STD -M $m0_iso $pm_iso $data_mt jmhat_mt_iso_on.dat $ctrl_fwd
#
#mpirun -np $procs $hostfile ../Mod3DMT_STD -M $m0_vti $pm_vti $data_mt jmhat_mt_vti_on.dat $ctrl_fwd
#
################################################################
# JMULT_T
################################################################
#
mpirun -np $procs $hostfile ../Mod3DMT_STD -T $m0_iso $data_mt dsigma_mt_iso_on.rho $ctrl_fwd
#
#mpirun -np $procs $hostfile ../Mod3DMT_STD -T $m0_vti $data_mt dsigma_mt_vti_on.rho $ctrl_fwd
#
################################################################
# NLCG
################################################################
#
mpirun -np $procs $hostfile ../Mod3DMT_STD -i -m $m0_iso $data_mt $ctrl_nlcg $ctrl_fwd
#
#mpirun -np $procs $hostfile ../Mod3DMT_STD -I NLCG $m0_vti $data_mt $ctrl_nlcg $ctrl_fwd
#
