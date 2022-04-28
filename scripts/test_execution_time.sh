#
# COMPILE MODEM-OO SERIAL
bash Configure.3D_MT.OSU.GFortran MakefileSerial SERIAL test_FWD.f90
make -f MakefileSerial clean
make -f MakefileSerial
#
# MODEM-OO SMALLEST EXAMPLE
./test_FWD -f --data ../inputs/esol/de.dat --model ../inputs/esol/pr.ws --control ../inputs/Others/first_control_file.txt -pd pred_data_serial_smallest.dat
#
# MODEM-OO MEDIUM EXAMPLE PARALLEL
./test_FWD -f --data ../inputs/1st_Example/rFile_Data_MT_LP.dat --model ../inputs/1st_Example/rFile_Model_trimed_lower_Boundary_50km --control ../inputs/Others/first_control_file.txt -pd pred_data_lp_trim.dat
#
# MODEM-OO BIGGEST EXAMPLE PARALLEL
./test_FWD -f --data ../inputs/1st_Example/rFile_Data_MT_TIP --model ../inputs/1st_Example/rFile_Model --control ../inputs/Others/first_control_file.txt -pd pred_data_tipper_biggest.dat
#
#
# COMPILE MODEM-OO MPI
bash Configure.3D_MT.OSU.GFortran MakefileMPI MPI test_FWD_MPI.f90
make -f MakefileMPI clean
make -f MakefileMPI
#
#
# MODEM-OO MEDIUM EXAMPLE SERIAL
mpirun -np 9 ./test_FWD_MPI -f --data ../inputs/1st_Example/rFile_Data_MT_LP.dat --model ../inputs/1st_Example/rFile_Model_trimed_lower_Boundary_50km --control ../inputs/Others/first_control_file.txt -pd pred_data_mpi_tipper.dat
#
# MODEM-OO MEDIUM EXAMPLE PARALLEL
./test_FWD -f --data ../inputs/1st_Example/rFile_Data_MT_LP.dat --model ../inputs/1st_Example/rFile_Model_trimed_lower_Boundary_50km --control ../inputs/Others/first_control_file.txt -pd pred_data_mpi_tipper.dat

#
# COMPILE MODEM-ON
make -f Makefile_SP2_EM1D_CSEM
mv Mod3DMT Mod3DMT_SP2
mpirun -np 33 ./Mod3DMT_SP2 -F ../../modem-oo/inputs/1st_Example/rFile_Model ../../modem-oo/inputs/1st_Example/rFile_Data_MT_LP.dat predicted_data_working.dat wFile_Esol ../inputs/modem_on_control_file.txt
mpirun -np 33 ./Mod3DMT_SP2 -F ../../modem-oo/inputs/1st_Example/rFile_Model_trimed_lower_Boundary_50km ../../modem-oo/inputs/1st_Example/rFile_Data_MT_LP.dat predicted_data_working_trimed.dat wFile_Esol ../inputs/modem_on_control_file.txt
