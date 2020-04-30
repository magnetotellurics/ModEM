#!/bin/bash
#
# ARGUMENTS: NONE
#
output_folder=pipeline_$(date "+%Y_%m_%d_%H_%M_%S")
#
echo "####### START $output_folder #######"
#
mkdir -p bin
#
mkdir -p outputs/temp
#
echo "####### SUMMARY FOR $output_folder #######\n\n" >> outputs/temp/summary.txt
echo "#" >> outputs/temp/summary.txt
echo "#" >> outputs/temp/summary.txt
#
#
echo "####### START build_modem_on #######"
#
#
echo "build_modem_on = build_modem_on.txt" >> outputs/temp/summary.txt
#
bash scripts/build_modem_on.sh | tee outputs/temp/build_modem_on.txt
#
echo "#" >> outputs/temp/summary.txt
#
#
echo "####### START update_doxygen #######"
#
mkdir -p docs/
#
echo "update_doxygen = update_doxygen.txt" >> outputs/temp/summary.txt
#
bash scripts/update_doxygen.sh | tee outputs/temp/update_doxygen.txt
#
echo "#" >> outputs/temp/summary.txt
#
#
echo "####### START test_forward_std_qmr #######"
#
echo "test_forward STD QMR MT = update_doxygen.txt" >> outputs/temp/summary.txt
#
bash scripts/test_forward.sh bin/Mod3DMT_STD inputs_workbench/rFile_userCtrl_MT-QMR $(nproc)
#
echo "#" >> outputs/temp/summary.txt
#
echo "test_forward STD QMR MT = update_doxygen.txt" >> outputs/temp/summary.txt
#
bash scripts/test_forward.sh bin/Mod3DMT_STD inputs_workbench/rFile_userCtrl_CSEM-QMR $(nproc)
#
echo "test_forward STD QMR CSEM = ????.txt" >> outputs/temp/summary.txt
#
#
echo "####### START test_forward_std_bicg #######"
#
echo "update_doxygen = update_doxygen.txt" >> outputs/temp/summary.txt
echo "#" >> outputs/temp/summary.txt
#
bash scripts/test_forward.sh bin/Mod3DMT_STD inputs_workbench/rFile_userCtrl_MT-BICG $(nproc)
#
echo "update_doxygen = update_doxygen.txt" >> outputs/temp/summary.txt
echo "#" >> outputs/temp/summary.txt
#
bash scripts/test_forward.sh bin/Mod3DMT_STD inputs_workbench/rFile_userCtrl_CSEM-BICG $(nproc)
#
#
echo "####### START test_forward_sp2_bicg #######"
#
echo "update_doxygen = update_doxygen.txt" >> outputs/temp/summary.txt
echo "#" >> outputs/temp/summary.txt
#
bash scripts/test_forward.sh bin/Mod3DMT_SP2 inputs_workbench/rFile_userCtrl_MT-BICG $(nproc)
#
echo "update_doxygen = update_doxygen.txt" >> outputs/temp/summary.txt
echo "#" >> outputs/temp/summary.txt
#
bash scripts/test_forward.sh bin/Mod3DMT_SP2 inputs_workbench/rFile_userCtrl_CSEM-BICG $(nproc)
#
#
mv outputs/temp outputs/$output_folder
