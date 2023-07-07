#!/bin/bash
#
# INPUTS../../modem_inputs/Benchmarking/rFile_Model_VTI
#
procs=22
m0_iso="../../modem_inputs/Benchmarking/rFile_Model"
pm_iso="../../modem_inputs/Benchmarking/rFile_Model"
m0_vti="../../modem_inputs/Benchmarking/rFile_Model_VTI"
pm_vti="../../modem_inputs/Benchmarking/rFile_Model_VTI"
data_mt="../../modem_inputs/Benchmarking/rFile_Data"
data_csem="../../modem_inputs/Benchmarking/Data_Ey.dat"
ctrl_dip1d="../../modem_inputs/Benchmarking/oo_fwd_control_dip1d"
ctrl_em1d="../../modem_inputs/Benchmarking/oo_fwd_control_em1d"
ctrl_nlcg="../../modem_inputs/Benchmarking/oo_inv_control_nlcg"
ctrl_dcg="../../modem_inputs/Benchmarking/oo_inv_control_dcg"
hostfile="--hostfile ../../hostfile_$procs.txt"
#hostfile=""
#
################################################################
# FWD
################################################################
#
mpirun -np $procs $hostfile ./ModEM_MPI -f -m $pm_iso -d $data_mt -pd pred_mt_iso_oo.dat -es esol_mt_iso_oo.bin -cf $ctrl_em1d
#
mpirun -np $procs $hostfile ./ModEM_MPI -f -m $pm_vti -d $data_mt -pd pred_mt_vti_oo.dat -es esol_mt_vti_oo.bin -cf $ctrl_em1d
#
#mpirun -np $procs $hostfile ./ModEM_MPI -f -m $pm_iso -d $data_csem -pd pred_dip1d_oo.dat -es esol_dip1d_oo.bin -cf $ctrl_dip1d
#
#mpirun -np $procs $hostfile ./ModEM_MPI -f -m $pm_iso -d $data_csem -pd pred_em1d_iso_oo.dat -es esol_em1d_iso_oo.bin -cf $ctrl_em1d
#
#mpirun -np $procs $hostfile ./ModEM_MPI -f -m $pm_vti -d $data_csem -pd pred_em1d_vti_oo.dat -es esol_em1d_vti_oo.bin -cf $ctrl_em1d
#
################################################################
# JMULT
################################################################
#
mpirun -np $procs $hostfile ./ModEM_MPI -j -m $m0_iso -pm $pm_iso -d $data_mt -jm jmhat_mt_iso_oo.dat -cf $ctrl_em1d
#
mpirun -np $procs $hostfile ./ModEM_MPI -j -m $m0_vti -pm $pm_vti -d $data_mt -jm jmhat_mt_vti_oo.dat -cf $ctrl_em1d
#
#mpirun -np $procs $hostfile ./ModEM_MPI -j -m $m0_iso -pm $pm_iso -d $data_csem -jm jmhat_dip1d_oo.dat -cf $ctrl_dip1d
#
#mpirun -np $procs $hostfile ./ModEM_MPI -j -m $m0_iso -pm $pm_iso -d $data_csem -jm jmhat_em1d_iso_oo.dat -cf $ctrl_em1d
#
#mpirun -np $procs $hostfile ./ModEM_MPI -j -m $m0_vti -pm $pm_vti -d $data_csem -jm jmhat_em1d_vti_oo.dat -cf $ctrl_em1d
#
################################################################
# JMULT_T
################################################################
#
mpirun -np $procs $hostfile ./ModEM_MPI -jt -m $m0_iso -d pred_mt_iso_oo.dat -dm dsigma_mt_iso_oo.rho -cf $ctrl_em1d
#
mpirun -np $procs $hostfile ./ModEM_MPI -jt -m $m0_vti -d pred_mt_vti_oo.dat -dm dsigma_mt_vti_oo.rho -cf $ctrl_em1d
#
#mpirun -np $procs $hostfile ./ModEM_MPI -jt -m $m0_iso -d pred_dip1d_oo.dat -dm dsigma_dip1d_oo.rho -cf $ctrl_dip1d
#
#mpirun -np $procs $hostfile ./ModEM_MPI -jt -m $m0_iso -d pred_em1d_iso_oo.dat -dm dsigma_em1d_iso_oo.rho -cf $ctrl_em1d
#
#mpirun -np $procs $hostfile ./ModEM_MPI -jt -m $m0_vti -d pred_em1d_vti_oo.dat -dm dsigma_em1d_vti_oo.rho -cf $ctrl_em1d
#
################################################################
# NLCG
################################################################
#
mpirun -np $procs $hostfile ./ModEM_MPI -i -m $m0_iso -d pred_mt_iso_oo.dat -ci $ctrl_nlcg -cf $ctrl_em1d -o mt_iso_oo_nlcg
#
mpirun -np $procs $hostfile ./ModEM_MPI -i -m $m0_vti -d pred_mt_vti_oo.dat -ci $ctrl_nlcg -cf $ctrl_em1d -o mt_vti_oo_nlcg
#
#mpirun -np $procs $hostfile ./ModEM_MPI -i -m $m0_iso -d pred_dip1d_oo.dat -ci $ctrl_nlcg -cf $ctrl_dip1d -o dip1d_oo_nlcg
#
#mpirun -np $procs $hostfile ./ModEM_MPI -i -m $m0_iso -d pred_em1d_iso_oo.dat -ci $ctrl_nlcg -cf $ctrl_em1d -o em1d_iso_oo_nlcg
#
#mpirun -np $procs $hostfile ./ModEM_MPI -i -m $m0_vti -d pred_em1d_vti_oo.dat -ci $ctrl_nlcg -cf $ctrl_em1d -o em1d_vti_oo_nlcg
#
################################################################
# DCG
################################################################
#
mpirun -np $procs $hostfile ./ModEM_MPI -i -m $m0_iso -d pred_mt_iso_oo.dat -ci $ctrl_dcg -cf $ctrl_em1d -o mt_iso_oo_dcg
#
mpirun -np $procs $hostfile ./ModEM_MPI -i -m $m0_vti -d pred_mt_vti_oo.dat -ci $ctrl_dcg -cf $ctrl_em1d -o mt_vti_oo_dcg
#
#./ModEM -i -m $m0_iso -d pred_dip1d_oo.dat -ci $ctrl_dcg -cf $ctrl_dip1d -o dip1d_oo_dcg
#
#./ModEM -i -m $m0_iso -d pred_em1d_iso_oo.dat -ci $ctrl_dcg -cf $ctrl_em1d -o em1d_iso_oo_dcg
#
#./ModEM -i -m $m0_vti -d pred_em1d_vti_oo.dat -ci $ctrl_dcg -cf $ctrl_em1d -o em1d_vti_oo_dcg
#
