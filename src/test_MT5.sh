#
mkdir OO_MT5_MF_QMR
#
cd OO_MT5_MF_QMR
#
mpirun -np 13 --hostfile ../../../hostfile_13.txt ../ModEM_MPI -i -m ../../../modem_inputs/MT_5BLOCKS/m0.rho -d ../../../modem_inputs/MT_5BLOCKS/d_05.dat -ci ../../../modem_inputs/MT_5BLOCKS/oo_inv_control_nlcg -cf ../control_mf_qmr -o inv_mt5_sp_bicg
#
mkdir ../OO_MT5_MF_BICG
#
cd ../OO_MT5_MF_BICG
#
mpirun -np 13 --hostfile ../../../hostfile_13.txt ../ModEM_MPI -i -m ../../../modem_inputs/MT_5BLOCKS/m0.rho -d ../../../modem_inputs/MT_5BLOCKS/d_05.dat -ci ../../../modem_inputs/MT_5BLOCKS/oo_inv_control_nlcg -cf ../control_mf_bicg -o inv_mt5_sp_bicg
#
mkdir ../OO_MT5_MF_BICG
#
cd ../OO_MT5_SP_QMR
#
mpirun -np 13 --hostfile ../../../hostfile_13.txt ../ModEM_MPI -i -m ../../../modem_inputs/MT_5BLOCKS/m0.rho -d ../../../modem_inputs/MT_5BLOCKS/d_05.dat -ci ../../../modem_inputs/MT_5BLOCKS/oo_inv_control_nlcg -cf ../control_sp qmr -o inv_mt5_sp_bicg
#
mkdir ../OO_MT5_MF_BICG
#
cd ../OO_MT5_SP_BICG
#
mpirun -np 13 --hostfile ../../../hostfile_13.txt ../ModEM_MPI -i -m ../../../modem_inputs/MT_5BLOCKS/m0.rho -d ../../../modem_inputs/MT_5BLOCKS/d_05.dat -ci ../../../modem_inputs/MT_5BLOCKS/oo_inv_control_nlcg -cf ../control_sp_bicg -o inv_mt5_sp_bicg
#
cd ..
#