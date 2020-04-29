#!/bin/bash
#
# ARGUMENTS: NONE
#
#
echo "START LOCAL PIPELINE SCRIPT"
#
mkdir -p bin
#
mkdir -p outputs/temp
#
bash scripts/build_modem_on.sh | tee outputs/temp/build_modem_on.txt
#
#mkdir -p docs/
#
#bash scripts/update_doxygen.sh | tee outputs/temp/update_doxygen.txt
#
bash scripts/test_plain_std.sh nproc | tee outputs/temp/test_plain_std.txt
#
bash scripts/test_plain_sp2.sh nproc | tee outputs/temp/test_plain_sp2.txt
#
bash scripts/test_forward.sh bin/Mod3DMT_STD inputs/rFile_Model inputs/rFile_Data nproc
#
bash scripts/test_forward.sh bin/Mod3DMT_SP2 inputs/rFile_Model inputs/rFile_Data nproc
#
bash scripts/test_forward_modelling.sh bin/Mod3DMT_STD inputs/rFile_userCtrl nproc
#
output_folder=output_$(date "+%Y_%m_%d_%H_%M_%S")
#
mv outputs/temp outputs/$output_folder
#
echo "FINISH LOCAL PIPELINE SCRIPT"