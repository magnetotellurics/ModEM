#
#
mkdir MF_SG_QMR
#
cd MF_SG_QMR
#
mpirun -np 13 ../ModEM_MPI -f -m ../../../modem_inputs/MR_test/m.ws -d ../../../modem_inputs/MR_test/d0.dat -ci ../../../modem_inputs/Benchmarking/oo_inv_control_nlcg -cf ../control_mf_qmr -o inv | tee -a output.txt
#
mkdir ../MF_SP_QMR
#
cd ../MF_SP_QMR
#
mpirun -np 13 ../ModEM_MPI -f -m ../../../modem_inputs/MR_test/m.ws -d ../../../modem_inputs/MR_test/d0.dat -ci ../../../modem_inputs/Benchmarking/oo_inv_control_nlcg -cf ../control_sp_qmr -o inv | tee -a output.txt
#
mkdir ../MF_SP2_QMR
#
cd ../MF_SP2_QMR
#
mpirun -np 13 ../ModEM_MPI -f -m ../../../modem_inputs/MR_test/m.ws -d ../../../modem_inputs/MR_test/d0.dat -ci ../../../modem_inputs/Benchmarking/oo_inv_control_nlcg -cf ../control_sp2_qmr -o inv | tee -a output.txt
#
mkdir ../MF_SG_BICG
#
cd ../MF_SG_BICG
#
mpirun -np 13 ../ModEM_MPI -f -m ../../../modem_inputs/MR_test/m.ws -d ../../../modem_inputs/MR_test/d0.dat -ci ../../../modem_inputs/Benchmarking/oo_inv_control_nlcg -cf ../control_mf_bicg -o inv | tee -a output.txt
#
mkdir ../MF_SP_BICG
#
cd ../MF_SP_BICG
#
mpirun -np 13 ../ModEM_MPI -f -m ../../../modem_inputs/MR_test/m.ws -d ../../../modem_inputs/MR_test/d0.dat -ci ../../../modem_inputs/Benchmarking/oo_inv_control_nlcg -cf ../control_sp_bicg -o inv | tee -a output.txt
#
mkdir ../MF_SP2_BICG
#
cd ../MF_SP2_BICG
#
mpirun -np 13 ../ModEM_MPI -f -m ../../../modem_inputs/MR_test/m.ws -d ../../../modem_inputs/MR_test/d0.dat -ci ../../../modem_inputs/Benchmarking/oo_inv_control_nlcg -cf ../control_sp2_bicg -o inv | tee -a output.txt
#