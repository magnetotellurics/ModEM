#!/bin/sh
#
#PBS -q workq 
#
app=$1
model=$2
data=$3
ctrl_inv=$4
ctrl_fwd=$5
cov=$6
out_dir=$7
#
mkdir -p ${out_dir}
#
echo "Starting INVERSE MODEM at $(date "+%d/%m/%Y - %H:%M:%S")" | tee -a ${out_dir}/std_output.txt
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
mpirun -np ${np} -hostfile ${PBS_NODEFILE} $app -I NLCG $model $data $ctrl_inv $ctrl_fwd $cov | tee -a ${out_dir}/std_output.txt 
#
cd
#
echo "Finish INVERSE MODEM" | tee -a  | tee -a ${out_dir}/std_output.txt
#
