#
# MODEM-OO BIG EXAMPLE SERIAL
bash Configure.3D_MT.OSU.GFortran MakefileSerial SERIAL test_FWD.f90
make -f MakefileSerial clean
make -f MakefileSerial
./test_FWD -f --data ../inputs/1st_Example/rFile_Data_MT_LP_TIP.dat --model ../inputs/1st_Example/rFile_Model --control ../inputs/Others/first_control_file.txt
mv predicted_data.dat predicted_data_oo
./test_FWD -f --data ../inputs/1st_Example/rFile_Data_MT_LP_TIP.dat --model ../inputs/1st_Example/rFile_Model_trimed_lower_Boundary_50km --control ../inputs/Others/first_control_file.txt
mv predicted_data.dat predicted_data_oo_trimed50
#
# MODEM-OO SMALL EXAMPLE SERIAL
./Configure.3D_MT.OSU.GFortran MakefileSerial SERIAL test_FWD.f90
make -f MakefileSerial clean
make -f MakefileSerial
./test_FWD -f --data ../inputs/esol/de.dat --model ../inputs/Others/simple_2_blocks.cpr --control ../inputs/Others/first_control_file.txt
mv predicted_data.dat predicted_data_00_2b
#
# MODEM-OO BIG EXAMPLE PARALLEL
bash Configure.3D_MT.OSU.GFortran MakefileMPI MPI test_FWD_MPI.f90
make -f MakefileMPI clean
make -f MakefileMPI
mpirun -np 9 ./test_FWD_MPI -f --data ../inputs/1st_Example/rFile_Data_MT --model ../inputs/1st_Example/rFile_Model --control ../inputs/first_control_file.txt
mv predicted_data.dat predicted_data_oo
mpirun -np 9 ./test_FWD_MPI -f --data ../inputs/1st_Example/rFile_Data_MT --model ../inputs/1st_Example/rFile_Model_trimed_lower_Boundary_50km --control ../inputs/first_control_file.txt
mv predicted_data.dat predicted_data_oo_trimed50
#
# MODEM-OO SMALL EXAMPLE PARALLEL
./Configure.3D_MT.OSU.GFortran MakefileMPI MPI test_FWD_MPI.f90
make -f MakefileMPI clean
make -f MakefileMPI
mpirun -np 5 ./test_FWD_MPI -f --data ../inputs/esol/de.dat --model ../inputs/Others/simple_2_blocks.cpr --control ../inputs/Others/first_control_file.txt
#
# MODEM-ON
make -f Makefile_SP2_EM1D_CSEM
mv Mod3DMT Mod3DMT_SP2
mpirun -np 33 ./Mod3DMT_SP2 -F ../../modem-oo/inputs/1st_Example/rFile_Model ../../modem-oo/inputs/1st_Example/rFile_Data_MT_LP.dat predicted_data_working.dat wFile_Esol ../inputs/modem_on_control_file.txt
mpirun -np 33 ./Mod3DMT_SP2 -F ../../modem-oo/inputs/1st_Example/rFile_Model_trimed_lower_Boundary_50km ../../modem-oo/inputs/1st_Example/rFile_Data_MT_LP.dat predicted_data_working_trimed.dat wFile_Esol ../inputs/modem_on_control_file.txt
