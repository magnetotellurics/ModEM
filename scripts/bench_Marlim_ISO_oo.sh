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
ctrl_fwd="../../../modem_inputs/V_4H/oo_fwd_control_em1d"
ctrl_nlcg="../../../modem_inputs/V_4H/oo_inv_control_nlcg"
hostfile="--hostfile ../../../hostfile_$procs.txt"
#hostfile=""
#
################################################################
# FWD
################################################################
#
mpirun -np $procs $hostfile ../ModEM_MPI -f -m $pm_iso -d $data -pd pred_mt_iso_oo.dat -es esol_mt_iso_oo.bin -cf $ctrl_fwd
#
#mpirun -np $procs $hostfile ../ModEM_MPI -f -m $pm_vti -d $data -pd pred_mt_vti_oo.dat -es esol_mt_vti_oo.bin -cf $ctrl_fwd
#
################################################################
# JMULT
################################################################
#
mpirun -np $procs $hostfile ../ModEM_MPI -j -m $m0_iso -pm $pm_iso -d $data -jm jmhat_mt_iso_oo.dat -cf $ctrl_fwd
#
#mpirun -np $procs $hostfile ../ModEM_MPI -j -m $m0_vti -pm $pm_vti -d $data -jm jmhat_mt_vti_oo.dat -cf $ctrl_fwd
#
################################################################
# JMULT_T
################################################################
#
mpirun -np $procs $hostfile ../ModEM_MPI -jt -m $m0_iso -d $data -dm dsigma_mt_iso_oo.rho -cf $ctrl_fwd
#
#mpirun -np $procs $hostfile ../ModEM_MPI -jt -m $m0_vti -d $data -dm dsigma_mt_vti_oo.rho -cf $ctrl_fwd
#
################################################################
# NLCG
################################################################
#
mpirun -np $procs $hostfile ../ModEM_MPI -i -m $m0_iso -d data -ci $ctrl_nlcg -cf $ctrl_fwd -o mt_iso_oo_nlcg -c $cov
#
#mpirun -np $procs $hostfile ../ModEM_MPI -i -m $m0_vti $data -ci $ctrl_nlcg -cf $ctrl_fwd -o mt_vti_oo_nlcg
#
