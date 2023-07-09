#!/bin/bash
#
# INPUTS
#
procs=22
m0_iso="../docs/vti_results/model_ISO_h0.mod"
pm_iso="../docs/vti_results/model_ISO_h10.mod"
m0_vti="../docs/vti_results/model_VTI_h0_v0.mod"
pm_vti="../docs/vti_results/model_VTI_h10_v10.mod"
data_mt="../docs/vti_results/Data_Z.dat"
ctrl_fwd="../docs/vti_results/oo_inv_control_nlcg"
ctrl_nlcg="../docs/vti_results/oo_inv_control_nlcg"
ctrl_dcg="../docs/vti_results/oo_inv_control_dcg"
hostfile="--hostfile ../../hostfile_$procs.txt"
#hostfile=""
#
################################################################
# FWD
################################################################
#
#mpirun -np $procs $hostfile ./ModEM_MPI -f -m $pm_iso -d $data_mt -pd pred_mt_iso_oo.dat -es esol_mt_iso_oo.bin -cf $ctrl_fwd
#
mpirun -np $procs $hostfile ./ModEM_MPI -f -m $pm_vti -d $data_mt -pd pred_mt_vti_oo.dat -es esol_mt_vti_oo.bin -cf $ctrl_fwd
#
################################################################
# JMULT
################################################################
#
#mpirun -np $procs $hostfile ./ModEM_MPI -j -m $m0_iso -pm $pm_iso -d $data_mt -jm jmhat_mt_iso_oo.dat -cf $ctrl_fwd
#
#mpirun -np $procs $hostfile ./ModEM_MPI -j -m $m0_vti -pm $pm_vti -d $data_mt -jm jmhat_mt_vti_oo.dat -cf $ctrl_fwd
#
################################################################
# JMULT_T
################################################################
#
#mpirun -np $procs $hostfile ./ModEM_MPI -jt -m $m0_iso -d pred_mt_iso_oo.dat -dm dsigma_mt_iso_oo.rho -cf $ctrl_fwd
#
#mpirun -np $procs $hostfile ./ModEM_MPI -jt -m $m0_vti -d pred_mt_vti_oo.dat -dm dsigma_mt_vti_oo.rho -cf $ctrl_fwd
#
################################################################
# NLCG
################################################################
#
#mpirun -np $procs $hostfile ./ModEM_MPI -i -m $m0_iso -d pred_mt_iso_oo.dat -ci $ctrl_nlcg -cf $ctrl_fwd -o mt_iso_oo_nlcg
#
mpirun -np $procs $hostfile ./ModEM_MPI -i -m $m0_vti -d pred_mt_vti_oo.dat -ci $ctrl_nlcg -cf $ctrl_fwd -o mt_vti_oo_nlcg
#
################################################################
# DCG
################################################################
#
#mpirun -np $procs $hostfile ./ModEM_MPI -i -m $m0_iso -d pred_mt_iso_oo.dat -ci $ctrl_dcg -cf $ctrl_fwd -o mt_iso_oo_dcg
#
#mpirun -np $procs $hostfile ./ModEM_MPI -i -m $m0_vti -d pred_mt_vti_oo.dat -ci $ctrl_dcg -cf $ctrl_fwd -o mt_vti_oo_dcg
#
