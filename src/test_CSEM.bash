#!/bin/bash
# #
# mpirun -np 55 --hostfile ../../hostfile_55.txt ./ModEM -f -m ../../modem_inputs/Small_model/Starting_model_Small_v_4_h.rho -d ../../modem_inputs/Small_model/Data_on_Small_Grid_20%Err.dat -cf control_mf_sg_qmr.txt -pd pred_mf_sg_qmr.dat | tee -a pred_mf_sg_qmr
# #
# mpirun -np 55 --hostfile ../../hostfile_55.txt ./ModEM -f -m ../../modem_inputs/Small_model/Starting_model_Small_v_4_h.rho -d ../../modem_inputs/Small_model/Data_on_Small_Grid_20%Err.dat -cf control_mf_sg_bicg.txt -pd pred_mf_sg_bicg.dat | tee -a pred_mf_sg_bicg
# #
# mpirun -np 55 --hostfile ../../hostfile_55.txt ./ModEM -f -m ../../modem_inputs/Small_model/Starting_model_Small_v_4_h.rho -d ../../modem_inputs/Small_model/Data_on_Small_Grid_20%Err.dat -cf control_sp_sg_qmr.txt -pd pred_sp_sg_qmr.dat | tee -a pred_sp_sg_qmr
# #
# mpirun -np 55 --hostfile ../../hostfile_55.txt ./ModEM -f -m ../../modem_inputs/Small_model/Starting_model_Small_v_4_h.rho -d ../../modem_inputs/Small_model/Data_on_Small_Grid_20%Err.dat -cf control_sp_sg_bicg.txt -pd pred_sp_sg_bicg.dat | tee -a pred_sp_sg_bicg
# #
# mpirun -np 55 --hostfile ../../hostfile_55.txt ./ModEM -f -m ../../modem_inputs/Small_model/Starting_model_Small_v_4_h.rho -d ../../modem_inputs/Small_model/Data_on_Small_Grid_20%Err.dat -cf control_sp2_sg_qmr.txt -pd pred_sp2_sg_qmr.dat | tee -a pred_sp2_sg_qmr
# #
# mpirun -np 55 --hostfile ../../hostfile_55.txt ./ModEM -f -m ../../modem_inputs/Small_model/Starting_model_Small_v_4_h.rho -d ../../modem_inputs/Small_model/Data_on_Small_Grid_20%Err.dat -cf control_sp2_sg_bicg.txt -pd pred_sp2_sg_bicg.dat | tee -a pred_sp2_sg_bicg
#
mpirun -np 55 --hostfile ../../hostfile_55.txt ./ModEM -i -m ../../modem_inputs/Small_model/Starting_model_Small_v_4_h.rho -d ../../modem_inputs/Small_model/Data_on_Small_Grid_20%Err.dat -ci inv_ctrl_template.txt -cf control_mf_sg_qmr.txt -o pred_mf_sg_qmr | tee -a inv_mf_sg_qmr
#
mpirun -np 55 --hostfile ../../hostfile_55.txt ./ModEM -i -m ../../modem_inputs/Small_model/Starting_model_Small_v_4_h.rho -d ../../modem_inputs/Small_model/Data_on_Small_Grid_20%Err.dat -ci inv_ctrl_template.txt -cf control_mf_sg_bicg.txt -o pred_mf_sg_bicg | tee -a inv_mf_sg_bicg
#
mpirun -np 55 --hostfile ../../hostfile_55.txt ./ModEM -i -m ../../modem_inputs/Small_model/Starting_model_Small_v_4_h.rho -d ../../modem_inputs/Small_model/Data_on_Small_Grid_20%Err.dat -ci inv_ctrl_template.txt -cf control_sp_sg_qmr.txt -o pred_sp_sg_qmr.dat | tee -a inv_sp_sg_qmr
#
mpirun -np 55 --hostfile ../../hostfile_55.txt ./ModEM -i -m ../../modem_inputs/Small_model/Starting_model_Small_v_4_h.rho -d ../../modem_inputs/Small_model/Data_on_Small_Grid_20%Err.dat -ci inv_ctrl_template.txt -cf control_sp_sg_bicg.txt -o pred_sp_sg_bicg.dat | tee -a inv_sp_sg_bicg
#
mpirun -np 55 --hostfile ../../hostfile_55.txt ./ModEM -i -m ../../modem_inputs/Small_model/Starting_model_Small_v_4_h.rho -d ../../modem_inputs/Small_model/Data_on_Small_Grid_20%Err.dat -ci inv_ctrl_template.txt -cf control_sp2_sg_qmr.txt -o pred_sp2_sg_qmr.dat | tee -a inv_sp2_sg_qmr
#
mpirun -np 55 --hostfile ../../hostfile_55.txt ./ModEM -i -m ../../modem_inputs/Small_model/Starting_model_Small_v_4_h.rho -d ../../modem_inputs/Small_model/Data_on_Small_Grid_20%Err.dat -ci inv_ctrl_template.txt -cf control_sp2_sg_bicg.txt -o pred_sp2_sg_bicg.dat | tee -a inv_sp2_sg_bicg
#
