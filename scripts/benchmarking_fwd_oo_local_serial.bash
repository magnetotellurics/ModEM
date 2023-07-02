#!/bin/bash
#
################################################################
# FWD
################################################################
#
./ModEM_SERIAL -f -m ../docs/vti_results/model_ISO_h10.mod -d ../docs/vti_results/Data_Z.dat -pd pred_mt_iso_oo_serial.dat -es esol_mt_iso_oo_serial.bin -cf ../docs/vti_results/oo_fwd_control_em1d
#
./ModEM_SERIAL -f -m ../docs/vti_results/model_VTI_h10_v10.mod -d ../docs/vti_results/Data_Z.dat -pd pred_mt_vti_oo_serial.dat -es esol_mt_vti_oo_serial.bin -cf ../docs/vti_results/oo_fwd_control_em1d
#
./ModEM_SERIAL -f -m ../docs/vti_results/model_ISO_h10.mod -d ../docs/vti_results/Data_Ey.dat -pd pred_dip1d_oo_serial.dat -es esol_dip1d_oo_serial.bin -cf ../docs/vti_results/oo_fwd_control_dip1d
#
./ModEM_SERIAL -f -m ../docs/vti_results/model_ISO_h10.mod -d ../docs/vti_results/Data_Ey.dat -pd pred_em1d_iso_oo_serial.dat -es esol_em1d_iso_oo_serial.bin -cf ../docs/vti_results/oo_fwd_control_em1d
#
./ModEM_SERIAL -f -m ../docs/vti_results/model_VTI_h10_v10.mod -d ../docs/vti_results/Data_Ey.dat -pd pred_em1d_vti_oo_serial.dat -es esol_em1d_vti_oo_serial.bin -cf ../docs/vti_results/oo_fwd_control_em1d
#
################################################################
# JMULT
################################################################
#
./ModEM_SERIAL -j -m ../docs/vti_results/model_ISO_h0.mod -pm ../docs/vti_results/model_ISO_h10.mod -d ../docs/vti_results/Data_Z.dat -jm jmhat_mt_iso_oo_serial.dat -cf ../docs/vti_results/oo_fwd_control_em1d
#
./ModEM_SERIAL -j -m ../docs/vti_results/model_VTI_h0_v0.mod -pm ../docs/vti_results/model_VTI_h10_v10.mod -d ../docs/vti_results/Data_Z.dat -jm jmhat_mt_vti_oo_serial.dat -cf ../docs/vti_results/oo_fwd_control_em1d
#
./ModEM_SERIAL -j -m ../docs/vti_results/model_ISO_h0.mod -pm ../docs/vti_results/model_ISO_h10.mod -d ../docs/vti_results/Data_Ey.dat -jm jmhat_dip1d_oo_serial.dat -cf ../docs/vti_results/oo_fwd_control_dip1d
#
./ModEM_SERIAL -j -m ../docs/vti_results/model_ISO_h0.mod -pm ../docs/vti_results/model_ISO_h10.mod -d ../docs/vti_results/Data_Ey.dat -jm jmhat_em1d_iso_oo_serial.dat -cf ../docs/vti_results/oo_fwd_control_em1d
#
./ModEM_SERIAL -j -m ../docs/vti_results/model_VTI_h0_v0.mod -pm ../docs/vti_results/model_VTI_h10_v10.mod -d ../docs/vti_results/Data_Ey.dat -jm jmhat_em1d_vti_oo_serial.dat -cf ../docs/vti_results/oo_fwd_control_em1d
#
################################################################
# JMULT_T
################################################################
#
./ModEM_SERIAL -jt -m ../docs/vti_results/model_ISO_h0.mod -d pred_mt_iso_oo_serial.dat -dm dsigma_mt_iso_oo_serial.rho -cf ../docs/vti_results/oo_fwd_control_em1d
#
./ModEM_SERIAL -jt -m ../docs/vti_results/model_VTI_h0_v0.mod -d pred_mt_vti_oo_serial.dat -dm dsigma_mt_vti_oo_serial.rho -cf ../docs/vti_results/oo_fwd_control_em1d
#
./ModEM_SERIAL -jt -m ../docs/vti_results/model_ISO_h0.mod -d pred_dip1d_oo_serial.dat -dm dsigma_dip1d_oo_serial.rho -cf ../docs/vti_results/oo_fwd_control_dip1d
#
./ModEM_SERIAL -jt -m ../docs/vti_results/model_ISO_h0.mod -d pred_em1d_iso_oo_serial.dat -dm dsigma_em1d_iso_oo_serial.rho -cf ../docs/vti_results/oo_fwd_control_em1d
#
./ModEM_SERIAL -jt -m ../docs/vti_results/model_VTI_h0_v0.mod -d pred_em1d_vti_oo_serial.dat -dm dsigma_em1d_vti_oo_serial.rho -cf ../docs/vti_results/oo_fwd_control_em1d
#
################################################################
# NLCG
################################################################
#
./ModEM_SERIAL -i -m ../docs/vti_results/model_ISO_h0.mod -d pred_mt_iso_oo_serial.dat -ci ../docs/vti_results/oo_inv_control_nlcg -cf ../docs/vti_results/oo_fwd_control_em1d -o mt_iso_oo_serial_nlcg
#
./ModEM_SERIAL -i -m ../docs/vti_results/model_VTI_h0_v0.mod -d pred_mt_vti_oo_serial.dat -ci ../docs/vti_results/oo_inv_control_nlcg -cf ../docs/vti_results/oo_fwd_control_em1d -o mt_vti_oo_serial_nlcg
#
./ModEM_SERIAL -i -m ../docs/vti_results/model_ISO_h0.mod -d pred_dip1d_oo_serial.dat -ci ../docs/vti_results/oo_inv_control_nlcg -cf ../docs/vti_results/oo_fwd_control_dip1d -o dip1d_oo_serial_nlcg
#
./ModEM_SERIAL -i -m ../docs/vti_results/model_ISO_h0.mod -d pred_em1d_iso_oo_serial.dat -ci ../docs/vti_results/oo_inv_control_nlcg -cf ../docs/vti_results/oo_fwd_control_em1d -o em1d_iso_oo_serial_nlcg
#
./ModEM_SERIAL -i -m ../docs/vti_results/model_VTI_h0_v0.mod -d pred_em1d_vti_oo_serial.dat -ci ../docs/vti_results/oo_inv_control_nlcg -cf ../docs/vti_results/oo_fwd_control_em1d -o em1d_vti_oo_serial_nlcg
#
################################################################
# DCG
################################################################
#
./ModEM_SERIAL -i -m ../docs/vti_results/model_ISO_h0.mod -d pred_mt_iso_oo_serial.dat -ci ../docs/vti_results/oo_inv_control_dcg -cf ../docs/vti_results/oo_fwd_control_em1d -o mt_iso_oo_serial_dcg
#
./ModEM_SERIAL -i -m ../docs/vti_results/model_VTI_h0_v0.mod -d pred_mt_vti_oo_serial.dat -ci ../docs/vti_results/oo_inv_control_dcg -cf ../docs/vti_results/oo_fwd_control_em1d -o mt_vti_oo_serial_dcg
#
./ModEM_SERIAL -i -m ../docs/vti_results/model_ISO_h0.mod -d pred_dip1d_oo_serial.dat -ci ../docs/vti_results/oo_inv_control_dcg -cf ../docs/vti_results/oo_fwd_control_dip1d -o dip1d_oo_serial_dcg
#
./ModEM_SERIAL -i -m ../docs/vti_results/model_ISO_h0.mod -d pred_em1d_iso_oo_serial.dat -ci ../docs/vti_results/oo_inv_control_dcg -cf ../docs/vti_results/oo_fwd_control_em1d -o em1d_iso_oo_serial_dcg
#
./ModEM_SERIAL -i -m ../docs/vti_results/model_VTI_h0_v0.mod -d pred_em1d_vti_oo_serial.dat -ci ../docs/vti_results/oo_inv_control_dcg -cf ../docs/vti_results/oo_fwd_control_em1d -o em1d_vti_oo_serial_dcg
#
