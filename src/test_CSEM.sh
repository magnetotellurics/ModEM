#
mkdir OO_CSEM20_MF_QMR
#
cd OO_CSEM20_MF_QMR
#
mpirun -np 55 --hostfile ../../../hostfile_55.txt ./ModEM_MPI -i -m ../../modem_inputs/CSEM_Small_model/Starting_model_Small_v_4_h.rho -d ../../modem_inputs/CSEM_Small_model/Data_on_Small_Grid_20%Err.dat -c ../../modem_inputs/CSEM_Small_model/Starting_model_small.cov -ci ../../modem_inputs/CSEM_Small_model/oo_inv_control_nlcg -cf control_mf_qmr -o inv_csem_mf_qmr
#
mkdir ../OO_CSEM20_MF_BICG
#
cd ../OO_CSEM20_MF_BICG
#
mpirun -np 55 --hostfile ../../../hostfile_55.txt ./ModEM_MPI -i -m ../../modem_inputs/CSEM_Small_model/Starting_model_Small_v_4_h.rho -d ../../modem_inputs/CSEM_Small_model/Data_on_Small_Grid_20%Err.dat -c ../../modem_inputs/CSEM_Small_model/Starting_model_small.cov -ci ../../modem_inputs/CSEM_Small_model/oo_inv_control_nlcg -cf control_mf_bicg -o inv_csem_mf_bicg
#
mkdir ../OO_CSEM20_MF_BICG
#
cd ../OO_CSEM20_SP_QMR
#
mpirun -np 55 --hostfile ../../../hostfile_55.txt ./ModEM_MPI -i -m ../../modem_inputs/CSEM_Small_model/Starting_model_Small_v_4_h.rho -d ../../modem_inputs/CSEM_Small_model/Data_on_Small_Grid_20%Err.dat -c ../../modem_inputs/CSEM_Small_model/Starting_model_small.cov -ci ../../modem_inputs/CSEM_Small_model/oo_inv_control_nlcg -cf control_sp_qmr -o inv_csem_sp_qmr
#
mkdir ../OO_CSEM20_MF_BICG
#
cd ../OO_CSEM20_SP_BICG
#
mpirun -np 55 --hostfile ../../../hostfile_55.txt ./ModEM_MPI -i -m ../../modem_inputs/CSEM_Small_model/Starting_model_Small_v_4_h.rho -d ../../modem_inputs/CSEM_Small_model/Data_on_Small_Grid_20%Err.dat -c ../../modem_inputs/CSEM_Small_model/Starting_model_small.cov -ci ../../modem_inputs/CSEM_Small_model/oo_inv_control_nlcg -cf control_sp_bicg -o inv_csem_sp_bicg
#
cd ..
#