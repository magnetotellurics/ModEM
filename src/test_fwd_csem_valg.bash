#!/bin/bash
#
# mf_sg_qmr
#
mpirun -np 55 ./ModEM -f -m ../../modem_inputs/Small_model/Starting_model_Small_v_4_h.rho -d ../../modem_inputs/Small_model/Data_on_Small_Grid_20%Err.dat -cf CTRL/fwd_control_mf_sg_qmr.txt -pd fwd_csem_mf_sg_qmr.dat | tee -a fwd_csem_mf_sg_qmr.txt
#
# mf_sg_bicg
#
mpirun -np 55 ./ModEM -f -m ../../modem_inputs/Small_model/Starting_model_Small_v_4_h.rho -d ../../modem_inputs/Small_model/Data_on_Small_Grid_20%Err.dat -cf CTRL/fwd_control_mf_sg_bicg.txt -pd fwd_csem_mf_sg_bicg.dat | tee -a fwd_csem_mf_sg_bicg.txt
#
# sp_sg_qmr
#
mpirun -np 55 ./ModEM -f -m ../../modem_inputs/Small_model/Starting_model_Small_v_4_h.rho -d ../../modem_inputs/Small_model/Data_on_Small_Grid_20%Err.dat -cf CTRL/fwd_control_sp_sg_qmr.txt -pd fwd_csem_sp_sg_qmr.dat | tee -a fwd_csem_sp_sg_qmr.txt
#
# sp_sg_bicg
#
mpirun -np 55 ./ModEM -f -m ../../modem_inputs/Small_model/Starting_model_Small_v_4_h.rho -d ../../modem_inputs/Small_model/Data_on_Small_Grid_20%Err.dat -cf CTRL/fwd_control_sp_sg_bicg.txt -pd fwd_csem_sp_sg_bicg.dat | tee -a fwd_csem_sp_sg_bicg.txt
#
# sp2_sg_qmr
#
mpirun -np 55 ./ModEM -f -m ../../modem_inputs/Small_model/Starting_model_Small_v_4_h.rho -d ../../modem_inputs/Small_model/Data_on_Small_Grid_20%Err.dat -cf CTRL/fwd_control_sp2_sg_qmr.txt -pd fwd_csem_sp2_sg_qmr.dat | tee -a fwd_csem_sp2_sg_qmr.txt
#
# sp2_sg_bicg
#
mpirun -np 55 ./ModEM -f -m ../../modem_inputs/Small_model/Starting_model_Small_v_4_h.rho -d ../../modem_inputs/Small_model/Data_on_Small_Grid_20%Err.dat -cf CTRL/fwd_control_sp2_sg_bicg.txt -pd fwd_csem_sp2_sg_bicg.dat | tee -a fwd_csem_sp2_sg_bicg.txt
#