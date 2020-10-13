#!/bin/bash
#
# ARGUMENTS: NONE
#
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
echo "build_modem_on = build_modem_on.txt" >> outputs/temp/summary.txt
START=$(date +%s%3N)
bash scripts/build_modem_on.sh | tee outputs/temp/build_modem_on.txt
END=$(date +%s%3N)
echo "Time Spent: $(( ( $END - $START ) / 1000 )) seconds" >> outputs/temp/summary.txt
echo "#" >> outputs/temp/summary.txt
#
#
echo "####### START update_doxygen #######"
#
mkdir -p docs/
#
echo "update_doxygen = update_doxygen.txt" >> outputs/temp/summary.txt
START=$(date +%s%3N)
bash scripts/update_doxygen.sh | tee outputs/temp/update_doxygen.txt
END=$(date +%s%3N)
echo "Time Spent: $(( ( $END - $START ) / 1000 )) seconds" >> outputs/temp/summary.txt
echo "#" >> outputs/temp/summary.txt
#
#
echo "####### START test_forward_std_qmr #######"
#
echo "test_forward STD QMR MT = ????.txt" >> outputs/temp/summary.txt
START=$(date +%s%3N)
bash scripts/test_forward.sh bin/Mod3DMT_STD inputs_workbench/rFile_userCtrl_MT-QMR $(nproc)
END=$(date +%s%3N)
echo "Time Spent: $(( ( $END - $START ) / 1000 )) seconds" >> outputs/temp/summary.txt
echo "#" >> outputs/temp/summary.txt
#
#
echo "test_forward STD QMR CSEM = ????.txt" >> outputs/temp/summary.txt
START=$(date +%s%3N)
bash scripts/test_forward.sh bin/Mod3DMT_STD inputs_workbench/rFile_userCtrl_CSEM-QMR $(nproc)
END=$(date +%s%3N)
echo "Time Spent: $(( ( $END - $START ) / 1000 )) seconds" >> outputs/temp/summary.txt
echo "#" >> outputs/temp/summary.txt
#
#
echo "####### START test_forward_std_bicg #######"
#
echo "test_forward STD MT BICG = ????.txt" >> outputs/temp/summary.txt
START=$(date +%s%3N)
bash scripts/test_forward.sh bin/Mod3DMT_STD inputs_workbench/rFile_userCtrl_MT-BICG $(nproc)
END=$(date +%s%3N)
echo "Time Spent: $(( ( $END - $START ) / 1000 )) seconds" >> outputs/temp/summary.txt
echo "#" >> outputs/temp/summary.txt
#
#
echo "test_forward STD CSEM BICG = ????.txt" >> outputs/temp/summary.txt
START=$(date +%s%3N)
bash scripts/test_forward.sh bin/Mod3DMT_STD inputs_workbench/rFile_userCtrl_CSEM-BICG $(nproc)
END=$(date +%s%3N)
echo "Time Spent: $(( ( $END - $START ) / 1000 )) seconds" >> outputs/temp/summary.txt
echo "#" >> outputs/temp/summary.txt
#
#
echo "####### START test_forward_sp2_bicg #######"#
#
echo "test_forward SP2 MT BICG = ????.txt" >> outputs/temp/summary.txt
START=$(date +%s%3N)
bash scripts/test_forward.sh bin/Mod3DMT_SP2 inputs_workbench/rFile_userCtrl_MT-BICG $(nproc)
END=$(date +%s%3N)
echo "Time Spent: $(( ( $END - $START ) / 1000 )) seconds" >> outputs/temp/summary.txt
echo "#" >> outputs/temp/summary.txt
#
#
echo "test_forward SP2 CSEM BICG = ????.txt" >> outputs/temp/summary.txt
START=$(date +%s%3N)
bash scripts/test_forward.sh bin/Mod3DMT_SP2 inputs_workbench/rFile_userCtrl_CSEM-BICG $(nproc)
END=$(date +%s%3N)
echo "Time Spent: $(( ( $END - $START ) / 1000 )) seconds" >> outputs/temp/summary.txt
echo "#" >> outputs/temp/summary.txt
#
#
mv outputs/temp outputs/$output_folder
