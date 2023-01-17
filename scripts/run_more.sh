#!/bin/bash
#
doxygen ../docs/doxygen_modem_oo_config
#
# MODEM-OO MT
#
#./Configure.modem.oo.GFortran MakefileSerial SERIAL TestSerial.f90
make -f MakefileSerial clean
make -f MakefileSerial
#
valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=csem_valgrind_fwd_serial.txt ./TestSerial -f -m ../inputs/1st_Example/rFile_Model_trim -d ../inputs/1st_Example/rFile_Data_CSEM -pd csem_fwd_pred_serial.txt -c ../docs/control_file_template
#
valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=csem_valgrind_adjt_serial.txt ./TestSerial -jt -m ../inputs/1st_Example/rFile_Model_trim -d ../inputs/1st_Example/rFile_Data_CSEM -pd csem_inv_pred_serial.txt -c ../docs/control_file_template
#
valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=mt4_valgrind_fwd_serial.txt ./TestSerial -f -m ../inputs/Benchmarking/rFile_Model -d ../inputs/Benchmarking/rFile_Data -pd mt4_fwd_pred_serial.txt -c ../docs/control_file_template
#
valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=mt4_valgrind_adjt_serial.txt ./TestSerial -jt -m ../inputs/Benchmarking/rFile_Model -d ../inputs/Benchmarking/rFile_Data -pd mt4_inv_pred_serial.txt -c ../docs/control_file_template
#

#./Configure.modem.oo.GFortran MakefileMPI MPI TestMPI.f90
make -f MakefileMPI clean
make -f MakefileMPI
#
#mpirun -np 5 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=mt4_valgrind_fwd_mpi.txt ./TestMPI -f -m ../inputs/Benchmarking/rFile_Model -d ../inputs/Benchmarking/rFile_Data -pd mt4_pred_mpi.txt -c ../docs/control_file_template
#
mpirun -np 5 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=mt4_valgrind_fwd_mpi.txt ./TestMPI -f -m ../inputs/Benchmarking/rFile_Model -d ../inputs/Benchmarking/rFile_Data -pd mt4_pred_mpi.txt -c ../docs/control_file_template
#
mkdir -p OO_MT
mv massif* OO_MT
#
# MODEM-ON MT
#
cd ../../modem-on/src
#
doxygen ../../modem-oo/docs/doxygen_modem_oo_config
#
mpirun -np 5 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=mt4_valgrind_fwd_mpi.txt ./Mod3DMT -F ../../modem-oo/inputs/Benchmarking/rFile_Model ../../modem-oo/inputs/Benchmarking/rFile_Data on_mt4_pred wfile ../../modem-oo/inputs/Others/modem_on_control_file.txt
#
mkdir -p ON_MT
mv massif* ON_MT
#
# MODEM-ON CSEM
#
mpirun -np 2 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=on_csem_valgrind_fwd_mpi.txt ./Mod3DMT -F ../../modem-oo/inputs/1st_Example/rFile_Model_trim ../../modem-oo/inputs/1st_Example/rFile_Data_CSEM on_csem_pred wfile ../../modem-oo/inputs/Others/modem_on_control_file.txt
#
mkdir -p ON_CSEM
mv massif* ON_CSEM
#
cd ../../modem-oo/src
#
mpirun -np 2 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=csem_valgrind_fwd_mpi.txt ./TestMPI -f -m ../inputs/1st_Example/rFile_Model_trim -d ../inputs/1st_Example/rFile_Data_CSEM -pd csem_mpi.txt -c ../docs/control_file_template
#
mkdir -p OO_CSEM
mv massif* OO_CSEM
#