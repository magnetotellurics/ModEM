#!/bin/bash
#
# INPUTS
#
procs=2
m0_iso="../docs/vti_results/model_ISO_h0.mod"
pm_iso="../docs/vti_results/model_ISO_h10.mod"
m0_vti="../docs/vti_results/model_VTI_h0_v0.mod"
pm_vti="../docs/vti_results/model_VTI_h10_v10.mod"
data_mt="../docs/vti_results/Data_Z_1tx_1rx.dat"
data_csem="../docs/vti_results/Data_Ey_1tx_1rx.dat"
ctrl_dip1d="../docs/vti_results/oo_fwd_control_dip1d"
ctrl_em1d="../docs/vti_results/oo_fwd_control_em1d"
ctrl_nlcg="../docs/vti_results/oo_inv_control_nlcg"
ctrl_dcg="../docs/vti_results/oo_inv_control_dcg"
#hostfile="--hostfile ../../hostfile_$procs.txt"
hostfile=""
#
################################################################
# FWD
################################################################
#
mpirun -np $procs $hostfile valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=pred_mt_iso_oo.vg ./ModEM_MPI -f -m $pm_iso -d $data_mt -pd pred_mt_iso_oo.dat -es esol_mt_iso_oo.bin -cf $ctrl_em1d
#
mpirun -np $procs $hostfile valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=pred_mt_iso_oo.vg ./ModEM_MPI -f -m $pm_vti -d $data_mt -pd pred_mt_vti_oo.dat -es esol_mt_vti_oo.bin -cf $ctrl_em1d
#
mpirun -np $procs $hostfile valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=pred_mt_iso_oo.vg ./ModEM_MPI -f -m $pm_iso -d $data_csem -pd pred_dip1d_oo.dat -es esol_dip1d_oo.bin -cf $ctrl_dip1d
#
mpirun -np $procs $hostfile valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=pred_mt_iso_oo.vg ./ModEM_MPI -f -m $pm_iso -d $data_csem -pd pred_em1d_iso_oo.dat -es esol_em1d_iso_oo.bin -cf $ctrl_em1d
#
mpirun -np $procs $hostfile valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=pred_mt_iso_oo.vg ./ModEM_MPI -f -m $pm_vti -d $data_csem -pd pred_em1d_vti_oo.dat -es esol_em1d_vti_oo.bin -cf $ctrl_em1d
#
################################################################
# JMULT
################################################################
#
mpirun -np $procs $hostfile valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=pred_mt_iso_oo.vg ./ModEM_MPI -j -m $m0_iso -pm $pm_iso -d $data_mt -jm jmhat_mt_iso_oo.dat -cf $ctrl_em1d
#
mpirun -np $procs $hostfile valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=pred_mt_iso_oo.vg ./ModEM_MPI -j -m $m0_vti -pm $pm_vti -d $data_mt -jm jmhat_mt_vti_oo.dat -cf $ctrl_em1d
#
mpirun -np $procs $hostfile valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=pred_mt_iso_oo.vg ./ModEM_MPI -j -m $m0_iso -pm $pm_iso -d $data_csem -jm jmhat_dip1d_oo.dat -cf $ctrl_dip1d
#
mpirun -np $procs $hostfile valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=pred_mt_iso_oo.vg ./ModEM_MPI -j -m $m0_iso -pm $pm_iso -d $data_csem -jm jmhat_em1d_iso_oo.dat -cf $ctrl_em1d
#
mpirun -np $procs $hostfile valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=pred_mt_iso_oo.vg ./ModEM_MPI -j -m $m0_vti -pm $pm_vti -d $data_csem -jm jmhat_em1d_vti_oo.dat -cf $ctrl_em1d
#
################################################################
# JMULT_T
################################################################
#
mpirun -np $procs $hostfile valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=pred_mt_iso_oo.vg ./ModEM_MPI -jt -m $m0_iso -d pred_mt_iso_oo.dat -dm dsigma_mt_iso_oo.rho -cf $ctrl_em1d
#
mpirun -np $procs $hostfile valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=pred_mt_iso_oo.vg ./ModEM_MPI -jt -m $m0_vti -d pred_mt_vti_oo.dat -dm dsigma_mt_vti_oo.rho -cf $ctrl_em1d
#
mpirun -np $procs $hostfile valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=pred_mt_iso_oo.vg ./ModEM_MPI -jt -m $m0_iso -d pred_dip1d_oo.dat -dm dsigma_dip1d_oo.rho -cf $ctrl_dip1d
#
mpirun -np $procs $hostfile valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=pred_mt_iso_oo.vg ./ModEM_MPI -jt -m $m0_iso -d pred_em1d_iso_oo.dat -dm dsigma_em1d_iso_oo.rho -cf $ctrl_em1d
#
mpirun -np $procs $hostfile valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=pred_mt_iso_oo.vg ./ModEM_MPI -jt -m $m0_vti -d pred_em1d_vti_oo.dat -dm dsigma_em1d_vti_oo.rho -cf $ctrl_em1d
#
################################################################
# NLCG
################################################################
#
#mpirun -np $procs $hostfile valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=pred_mt_iso_oo.vg ./ModEM_MPI -i -m $m0_iso -d pred_mt_iso_oo.dat -ci $ctrl_nlcg -cf $ctrl_em1d -o mt_iso_oo_nlcg
#
#mpirun -np $procs $hostfile valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=pred_mt_iso_oo.vg ./ModEM_MPI -i -m $m0_vti -d pred_mt_vti_oo.dat -ci $ctrl_nlcg -cf $ctrl_em1d -o mt_vti_oo_nlcg
#
#mpirun -np $procs $hostfile valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=pred_mt_iso_oo.vg ./ModEM_MPI -i -m $m0_iso -d pred_dip1d_oo.dat -ci $ctrl_nlcg -cf $ctrl_dip1d -o dip1d_oo_nlcg
#
#mpirun -np $procs $hostfile valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=pred_mt_iso_oo.vg ./ModEM_MPI -i -m $m0_iso -d pred_em1d_iso_oo.dat -ci $ctrl_nlcg -cf $ctrl_em1d -o em1d_iso_oo_nlcg
#
#mpirun -np $procs $hostfile valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=pred_mt_iso_oo.vg ./ModEM_MPI -i -m $m0_vti -d pred_em1d_vti_oo.dat -ci $ctrl_nlcg -cf $ctrl_em1d -o em1d_vti_oo_nlcg
#
################################################################
# DCG
################################################################
#
#./ModEM -i -m $m0_iso -d pred_mt_iso_oo.dat -ci $ctrl_dcg -cf $ctrl_em1d -o mt_iso_oo_dcg
#
#./ModEM -i -m $m0_vti -d pred_mt_vti_oo.dat -ci $ctrl_dcg -cf $ctrl_em1d -o mt_vti_oo_dcg
#
#./ModEM -i -m $m0_iso -d pred_dip1d_oo.dat -ci $ctrl_dcg -cf $ctrl_dip1d -o dip1d_oo_dcg
#
#./ModEM -i -m $m0_iso -d pred_em1d_iso_oo.dat -ci $ctrl_dcg -cf $ctrl_em1d -o em1d_iso_oo_dcg
#
#./ModEM -i -m $m0_vti -d pred_em1d_vti_oo.dat -ci $ctrl_dcg -cf $ctrl_em1d -o em1d_vti_oo_dcg
#
