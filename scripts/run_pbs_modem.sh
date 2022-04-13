#!/bin/bash
#
app=$1
user_ctrl=$2
out_dir=$3
#
work_dir=$(pwd)
#
# GET JOB NAME FROM USER CTRL
grep -hnr "job" ${user_ctrl} > job_line.txt
job=$(sed '1q;d' job_line.txt | awk '{ print $2 }')
rm job_line.txt
#
# GET DATA FILE NAME FROM USER CTRL
grep -hnr "rFile_Data" ${user_ctrl} > data_line.txt
data=$(sed '1q;d' data_line.txt | awk '{ print $2 }')
rm data_line.txt
data=$(echo $data | cut -d "'" -f 2)
#
# AUTOMATED CALCULATION OF NODES|PPN
cd scripts/CalcProcs
g++ *.cpp -O3 -o CalcProcs
node_param=$(./CalcProcs $data)
rm CalcProcs
#
cd ${work_dir}
#
#echo "scripts/run_modem.sh ${app} ${user_ctrl} ${out_dir} | qsub -l ${node_param} -N $(hostname)_${job}_$(date "+%m%d%H%M") -o ${out_dir}/pbs_output.txt -e ${out_dir}/pbs_error.txt"
#
echo "scripts/run_modem.sh ${app} ${user_ctrl} ${out_dir}" | qsub -l ${node_param} -N $(hostname)_${job}_$(date "+%m%d%H%M") -o ${out_dir}/pbs_output.txt -e ${out_dir}/pbs_error.txt
#
