#!/bin/sh
#
#PBS -q workq 
#
app=$1
user_ctrl=$2
out_dir=$3
#
mkdir -p ${out_dir}
#
echo "Starting FWD MODEM at $(date "+%d/%m/%Y - %H:%M:%S")" | tee -a ${out_dir}/std_output.txt
#
np=$(wc -l < ${PBS_NODEFILE})
echo "${np} process" | tee -a ${out_dir}/std_output.txt
#
echo "PBS_NODEFILE : ${PBS_NODEFILE}" | tee -a ${out_dir}/std_output.txt
cat ${PBS_NODEFILE} | tee -a ${out_dir}/std_output.txt
#
PBS_WORKDIR=${out_dir}
#
cd ${out_dir}
#
mpirun -np ${np} -hostfile ${PBS_NODEFILE} $app -W ${user_ctrl} | tee -a ${out_dir}/std_output.txt 
#
cd
#
echo "Finish FWD MODEM" | tee -a  | tee -a ${out_dir}/std_output.txt
#
