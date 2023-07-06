#!/bin/bash
#
# INPUTS
#
m0_iso="../docs/vti_results/model_ISO_h0.mod"
pm_iso="../docs/vti_results/model_ISO_h10.mod"
m0_vti="../docs/vti_results/model_VTI_h0_v0.mod"
pm_vti="../docs/vti_results/model_VTI_h10_v10.mod"
data_mt="../docs/vti_results/Data_Z.dat"
data_csem="../docs/vti_results/Data_Ey.dat"
ctrl_dip1d="../docs/vti_results/oo_fwd_control_dip1d"
ctrl_em1d="../docs/vti_results/oo_fwd_control_em1d"
ctrl_nlcg="../docs/vti_results/oo_inv_control_nlcg"
ctrl_dcg="../docs/vti_results/oo_inv_control_dcg"
#
################################################################
# FWD
################################################################
#
./ModEM_SERIAL -f -m $pm_iso -d $data_mt -pd pred_mt_iso_oo_serial.dat -es esol_mt_iso_oo_serial.bin -cf $ctrl_em1d
#
./ModEM_SERIAL -f -m $pm_vti -d $data_mt -pd pred_mt_vti_oo_serial.dat -es esol_mt_vti_oo_serial.bin -cf $ctrl_em1d
#
./ModEM_SERIAL -f -m $pm_iso -d $data_csem -pd pred_dip1d_oo_serial.dat -es esol_dip1d_oo_serial.bin -cf $ctrl_dip1d
#
./ModEM_SERIAL -f -m $pm_iso -d $data_csem -pd pred_em1d_iso_oo_serial.dat -es esol_em1d_iso_oo_serial.bin -cf $ctrl_em1d
#
./ModEM_SERIAL -f -m $pm_vti -d $data_csem -pd pred_em1d_vti_oo_serial.dat -es esol_em1d_vti_oo_serial.bin -cf $ctrl_em1d
#
################################################################
# JMULT
################################################################
#
./ModEM_SERIAL -j -m $m0_iso -pm $pm_iso -d $data_mt -jm jmhat_mt_iso_oo_serial.dat -cf $ctrl_em1d
#
./ModEM_SERIAL -j -m $m0_vti -pm $pm_vti -d $data_mt -jm jmhat_mt_vti_oo_serial.dat -cf $ctrl_em1d
#
./ModEM_SERIAL -j -m $m0_iso -pm $pm_iso -d $data_csem -jm jmhat_dip1d_oo_serial.dat -cf $ctrl_dip1d
#
./ModEM_SERIAL -j -m $m0_iso -pm $pm_iso -d $data_csem -jm jmhat_em1d_iso_oo_serial.dat -cf $ctrl_em1d
#
./ModEM_SERIAL -j -m $m0_vti -pm $pm_vti -d $data_csem -jm jmhat_em1d_vti_oo_serial.dat -cf $ctrl_em1d
#
################################################################
# JMULT_T
################################################################
#
./ModEM_SERIAL -jt -m $m0_iso -d pred_mt_iso_oo_serial.dat -dm dsigma_mt_iso_oo_serial.rho -cf $ctrl_em1d
#
./ModEM_SERIAL -jt -m $m0_vti -d pred_mt_vti_oo_serial.dat -dm dsigma_mt_vti_oo_serial.rho -cf $ctrl_em1d
#
./ModEM_SERIAL -jt -m $m0_iso -d pred_dip1d_oo_serial.dat -dm dsigma_dip1d_oo_serial.rho -cf $ctrl_dip1d
#
./ModEM_SERIAL -jt -m $m0_iso -d pred_em1d_iso_oo_serial.dat -dm dsigma_em1d_iso_oo_serial.rho -cf $ctrl_em1d
#
./ModEM_SERIAL -jt -m $m0_vti -d pred_em1d_vti_oo_serial.dat -dm dsigma_em1d_vti_oo_serial.rho -cf $ctrl_em1d
#
################################################################
# NLCG
################################################################
#
./ModEM_SERIAL -i -m $m0_iso -d pred_mt_iso_oo_serial.dat -ci $ctrl_nlcg -cf $ctrl_em1d -o mt_iso_oo_serial_nlcg
#
./ModEM_SERIAL -i -m $m0_vti -d pred_mt_vti_oo_serial.dat -ci $ctrl_nlcg -cf $ctrl_em1d -o mt_vti_oo_serial_nlcg
#
./ModEM_SERIAL -i -m $m0_iso -d pred_dip1d_oo_serial.dat -ci $ctrl_nlcg -cf $ctrl_dip1d -o dip1d_oo_serial_nlcg
#
./ModEM_SERIAL -i -m $m0_iso -d pred_em1d_iso_oo_serial.dat -ci $ctrl_nlcg -cf $ctrl_em1d -o em1d_iso_oo_serial_nlcg
#
./ModEM_SERIAL -i -m $m0_vti -d pred_em1d_vti_oo_serial.dat -ci $ctrl_nlcg -cf $ctrl_em1d -o em1d_vti_oo_serial_nlcg
#
################################################################
# DCG
################################################################
#
#./ModEM_SERIAL -i -m $m0_iso -d pred_mt_iso_oo_serial.dat -ci $ctrl_dcg -cf $ctrl_em1d -o mt_iso_oo_serial_dcg
#
#./ModEM_SERIAL -i -m $m0_vti -d pred_mt_vti_oo_serial.dat -ci $ctrl_dcg -cf $ctrl_em1d -o mt_vti_oo_serial_dcg
#
#./ModEM_SERIAL -i -m $m0_iso -d pred_dip1d_oo_serial.dat -ci $ctrl_dcg -cf $ctrl_dip1d -o dip1d_oo_serial_dcg
#
#./ModEM_SERIAL -i -m $m0_iso -d pred_em1d_iso_oo_serial.dat -ci $ctrl_dcg -cf $ctrl_em1d -o em1d_iso_oo_serial_dcg
#
#./ModEM_SERIAL -i -m $m0_vti -d pred_em1d_vti_oo_serial.dat -ci $ctrl_dcg -cf $ctrl_em1d -o em1d_vti_oo_serial_dcg
#
