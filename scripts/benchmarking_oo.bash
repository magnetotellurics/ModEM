#!/bin/bash
#
nprocs=5
name='BLOCK2'
model='../../modem_inputs/Benchmarking/rFile_Model'
pmodel='../../modem_inputs/Benchmarking/rFile_Model'
data='../../modem_inputs/Benchmarking/rFile_Data'
#
mkdir ${name}_OO
mkdir ${name}_OO/SERIAL
mkdir ${name}_OO/MPI
#
# FWD
#
./ModEM_SERIAL -f -m $model -d $data -cf fwd_ctrl_template.txt -pd ${name}_pred.dat
mv ${name}_pred.dat ${name}_OO/SERIAL/
#
mpirun -np $nprocs ./ModEM_MPI -f -m $model -d $data -cf fwd_ctrl_template.txt -pd ${name}_pred.dat
mv ${name}_pred.dat ${name}_OO/MPI/
#
# JMULT
#
./ModEM_SERIAL -j -m $model -pm $pmodel -d $data -cf fwd_ctrl_template.txt -jm ${name}_jmhat.dat
mv ${name}_jmhat.dat ${name}_OO/SERIAL/
#
mpirun -np $nprocs ./ModEM_MPI -j -m $model -pm $pmodel -d $data -cf fwd_ctrl_template.txt -jm ${name}_jmhat.dat
mv ${name}_jmhat.dat ${name}_OO/MPI/
#
# JMULT_T
#
./ModEM_SERIAL -jt -m $model -d $data -cf fwd_ctrl_template.txt -dm ${name}_sigma.rho
mv ${name}_sigma.rho ${name}_OO/SERIAL/
#
mpirun -np $nprocs ./ModEM_MPI -jt -m $model -d $data -cf fwd_ctrl_template.txt -dm ${name}_sigma.rho
mv ${name}_sigma.rho ${name}_OO/MPI/
#
# INV
#
./ModEM_SERIAL -i -m $model -d $data -cf fwd_ctrl_template.txt -ci inv_ctrl_template.txt -o ${name}_OO/SERIAL/INV
mv fort.198* ${name}_OO/SERIAL/INV/
#
mpirun -np $nprocs ./ModEM_MPI -i -m $model -d $data -cf fwd_ctrl_template.txt -ci inv_ctrl_template.txt -o ${name}_OO/MPI/INV
mv fort.198* ${name}_OO/MPI/INV/
#