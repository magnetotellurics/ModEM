#!/bin/bash
#
app=$1
model=$2
data=$3
em_soln=$4
fwd_ctrl=$5
out_dir=$6
#
work_dir=$(pwd)
#
# AUTOMATED CALCULATION OF NODES|PPN
cd scripts/CalcProcs
g++ *.cpp -O3 -o CalcProcs
node_param=$(./CalcProcs $data)
rm CalcProcs
#
echo ${node_param}
#
cd ${work_dir}
#
echo "scripts/pbs/fwd_modem.sh $app $model $data $em_soln $fwd_ctrl $out_dir" | qsub -l ${node_param} -N FWD_$(date "+%m%d%H%M") -o ${out_dir}/pbs_output.txt -e ${out_dir}/pbs_error.txt
#
