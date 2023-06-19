#!/bin/bash
#
nprocs=9
name='BLOCK2'
model='../../modem_inputs/Benchmarking/rFile_Model'
pmodel='../../modem_inputs/Benchmarking/rFile_Model'
data='../../modem_inputs/Benchmarking/rFile_Data'
#
mkdir ${name}_ON
mkdir ${name}_ON/FWD
mkdir ${name}_ON/JMULT
mkdir ${name}_ON/JMULT_T
mkdir ${name}_ON/INV
#
# FWD
#
mpirun -np $nprocs ./Mod3DMT_STD -F $model $data ${name}_pred.dat esoln.bin control.fwd
mv Nodes_Status_* ${name}_ON/FWD/
mv QMR_* ${name}_ON/FWD/
mv ${name}_pred.dat ${name}_ON/FWD/
mv esoln ${name}_ON/FWD/
#
# JMULT
#
mpirun -np $nprocs ./Mod3DMT_STD -M $model $pmodel $data ${name}_jmhat.dat control.fwd
mv Nodes_Status_* ${name}_ON/JMULT/
mv QMR_* ${name}_ON/JMULT/
mv ${name}_jmhat.dat ${name}_ON/JMULT/
#
# JMULT_T
#
mpirun -np $nprocs ./Mod3DMT_STD -T $model $data ${name}_dsigma.rho control.fwd
mkdir JMULT_T
mv Nodes_Status_* ${name}_ON/JMULT_T/
mv QMR_* JMULT_T/
mv ${name}_dsigma.rho ${name}_ON/JMULT_T/
#
# INV
#
mpirun -np $nprocs ./Mod3DMT_STD -I NLCG $model $data control.inv control.fwd
mv ${name}* ${name}_ON/INV/
mv fort.198* ${name}_ON/INV/
mv Nodes_Status_* ${name}_ON/INV/
mv QMR_* ${name}_ON/INV/
#