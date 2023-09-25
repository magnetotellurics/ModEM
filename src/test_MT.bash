#!/bin/bash
# #
# mpirun -np 13 ./ModEM -f -m ../../modem_inputs/MR_test/m_fixed.ws -d ../../modem_inputs/MR_test/d0_fixed.dat -cf control_mf_sg_qmr.txt -pd pred_mf_sg_qmr.dat | tee -a pred_mf_sg_qmr
# #
# mpirun -np 13 ./ModEM -f -m ../../modem_inputs/MR_test/m_fixed.ws -d ../../modem_inputs/MR_test/d0_fixed.dat -cf control_mf_sg_bicg.txt -pd pred_mf_sg_bicg.dat | tee -a pred_mf_sg_bicg
# #
# mpirun -np 13 ./ModEM -f -m ../../modem_inputs/MR_test/m_fixed.ws -d ../../modem_inputs/MR_test/d0_fixed.dat -cf control_sp_sg_qmr.txt -pd pred_sp_sg_qmr.dat | tee -a pred_sp_sg_qmr
# #
# mpirun -np 13 ./ModEM -f -m ../../modem_inputs/MR_test/m_fixed.ws -d ../../modem_inputs/MR_test/d0_fixed.dat -cf control_sp_sg_bicg.txt -pd pred_sp_sg_bicg.dat | tee -a pred_sp_sg_bicg
# #
# mpirun -np 13 ./ModEM -f -m ../../modem_inputs/MR_test/m_fixed.ws -d ../../modem_inputs/MR_test/d0_fixed.dat -cf control_sp2_sg_qmr.txt -pd pred_sp2_sg_qmr.dat | tee -a pred_sp2_sg_qmr
# #
# mpirun -np 13 ./ModEM -f -m ../../modem_inputs/MR_test/m_fixed.ws -d ../../modem_inputs/MR_test/d0_fixed.dat -cf control_sp2_sg_bicg.txt -pd pred_sp2_sg_bicg.dat | tee -a pred_sp2_sg_bicg
#
mpirun -np 13 ./ModEM -i -m ../../modem_inputs/MR_test/m_fixed.ws -d ../../modem_inputs/MR_test/d0_fixed.dat -ci inv_ctrl_template.txt -cf control_mf_sg_qmr.txt -o pred_mf_sg_qmr | tee -a inv_mf_sg_qmr
#
mpirun -np 13 ./ModEM -i -m ../../modem_inputs/MR_test/m_fixed.ws -d ../../modem_inputs/MR_test/d0_fixed.dat -ci inv_ctrl_template.txt -cf control_mf_sg_bicg.txt -o pred_mf_sg_bicg | tee -a inv_mf_sg_bicg
#
mpirun -np 13 ./ModEM -i -m ../../modem_inputs/MR_test/m_fixed.ws -d ../../modem_inputs/MR_test/d0_fixed.dat -ci inv_ctrl_template.txt -cf control_sp_sg_qmr.txt -o pred_sp_sg_qmr.dat | tee -a inv_sp_sg_qmr
#
mpirun -np 13 ./ModEM -i -m ../../modem_inputs/MR_test/m_fixed.ws -d ../../modem_inputs/MR_test/d0_fixed.dat -ci inv_ctrl_template.txt -cf control_sp_sg_bicg.txt -o pred_sp_sg_bicg.dat | tee -a inv_sp_sg_bicg
#
mpirun -np 13 ./ModEM -i -m ../../modem_inputs/MR_test/m_fixed.ws -d ../../modem_inputs/MR_test/d0_fixed.dat -ci inv_ctrl_template.txt -cf control_sp2_sg_qmr.txt -o pred_sp2_sg_qmr.dat | tee -a inv_sp2_sg_qmr
#
mpirun -np 13 ./ModEM -i -m ../../modem_inputs/MR_test/m_fixed.ws -d ../../modem_inputs/MR_test/d0_fixed.dat -ci inv_ctrl_template.txt -cf control_sp2_sg_bicg.txt -o pred_sp2_sg_bicg.dat | tee -a inv_sp2_sg_bicg
#
