#!/bin/bash
#
echo "MAIN SCRIPT STARTS: $(date "+%Y/%m/%d - %H:%M:%S")" | tee -a summary.txt
#
T_START=$(date +%s%3N)
#
# FWD
#
mkdir TEST_CSEM_VALGRIND_FWD
	#
	cd TEST_CSEM_VALGRIND_FWD
		#
		valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=oo_valgrind_output_csem.txt ../test_FWD --forward --model ../../inputs/Naser_CSEM/rFile_Model  --data ../../inputs/Naser_CSEM/rFile_Data_1_fix --control ../../inputs/Others/first_control_file.txt | tee -a output.txt
		#
		T_END=$(date +%s%3N)
		echo "	> TEST_CSEM_VALGRIND_FWD: $(( ( $T_END - $T_START ) / 1000 )) seconds" | tee -a ../summary.txt
		T_START=$(date +%s%3N)
		#
	cd ..
	#
mkdir TEST_CSEM_MASSIF_FWD
	#
	cd TEST_CSEM_MASSIF_FWD
		#
		valgrind --tool=massif --time-unit=ms --max-snapshots=1000 ../test_FWD --forward --model ../../inputs/Naser_CSEM/rFile_Model  --data ../../inputs/Naser_CSEM/rFile_Data_1_fix --control ../../inputs/Others/first_control_file.txt | tee -a output.txt
		#
		bash ../../scripts/memory_brief.sh
		#
		T_END=$(date +%s%3N)
		echo "	> TEST_CSEM_MASSIF_FWD: $(( ( $T_END - $T_START ) / 1000 )) seconds" | tee -a ../summary.txt
		T_START=$(date +%s%3N)
		#
	cd ..
	#
echo "MAIN SCRIPT FINISHES: $(date "+%Y/%m/%d - %H:%M:%S")" | tee -a summary.txt
#
exit 0
#
# END OF SCRIPT
