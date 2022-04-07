#
#
./Configure.3D_MT.OSU.GFortran MakefileSerial SERIAL test_FWD.f90
make -f MakefileSerial clean
make -f MakefileSerial
./test_FWD -f --data ../inputs/esol/de.dat --model ../inputs/simple_2_blocks.cpr --control ../inputs/first_control_file.txt
#
#
./Configure.3D_MT.OSU.GFortran MakefileMPI MPI test_FWD_MPI.f90
make -f MakefileMPI clean
make -f MakefileMPI
mpirun -np 5 ./test_FWD_MPI -f --data ../inputs/esol/de.dat --model ../inputs/simple_2_blocks.cpr --control ../inputs/first_control_file.txt
