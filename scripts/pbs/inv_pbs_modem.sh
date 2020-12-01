#!/bin/bash
#
app=$1
model=$2
data=$3
ctrl_inv=$4
ctrl_fwd=$5
cov=$6
out_dir=$7
#
work_dir=$(pwd)
#
# AUTOMATED CALCULATION OF NODES|PPN
cd scripts/CalcProcs
g++ *.cpp -O3 -o CalcProcs
node_param=$(./CalcProcs $data)
rm CalcProcs
#
cd ${work_dir}
#
echo "scripts/inv_modem.sh $app $model $data $ctrl_inv $ctrl_fwd $cov $out_dir" | qsub -l ${node_param} -N INV_$(date "+%m%d%H%M") -o ${out_dir}/pbs_output.txt -e ${out_dir}/pbs_error.txt
#
