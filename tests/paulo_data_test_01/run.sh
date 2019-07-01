#!/bin/bash
#
# ARGUMENTS: 
#		$1 - FIRST DATA FULL PATH
#		$2 - SECOND DATA FULL PATH
#		$3 - OUTPUT DIR
#
data_1=$1
data_2=$2
output_dir=$3
#
emptyspace=""
# 
data_1_name=$data_1
data_1_name=${data_1_name/.dat/$emptyspace}
data_1_name=${data_1_name/*\//$emptyspace}
# 
data_2_name=$data_2
data_2_name=${data_2_name/.dat/$emptyspace}
data_2_name=${data_2_name/*\//$emptyspace}
#
echo "### START 'data_explorer' BETWEEN '$data_1_name' AND '$data_2_name' ###" >> std_out.txt
#
cd ../data_explorer/
#
#
g++ *.cpp -o data_explorer
#
#
./data_explorer ../$output_dir$data_1 $data_2 &>> ../$output_dir/cmp_results.txt
#
rm -rf data_explorer
#
cd ../$output_dir
#
echo "### FINISH 'data_explorer' ###" >> std_out.txt
#
#
exit 0
#
# END OF SCRIPT
