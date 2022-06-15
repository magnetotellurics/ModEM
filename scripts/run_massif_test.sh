#!/bin/bash
#
echo "MAIN SCRIPT STARTS: $(date "+%Y/%m/%d - %H:%M:%S")" | tee -a summary.txt
#
T_START=$(date +%s%3N)
#
# FWD
#
mkdir TEST_CSEM_MASSIF_FWD
	#
	cd TEST_CSEM_MASSIF_FWD
		#
		mpirun -np 2 valgrind --tool=massif --time-unit=ms --max-snapshots=1000 ../Mod3DMT -F ../../inputs/1st_Example/rFile_Model_trim ../../inputs/1st_Example/rFile_Data_CSEM predicted | tee -a output.txt
		#
		bash ../../scripts/memory_brief.sh
		#
		T_END=$(date +%s%3N)
		echo "	> TEST_CSEM_MASSIF_FWD: $(( ( $T_END - $T_START ) / 1000 )) seconds" | tee -a ../summary.txt
		T_START=$(date +%s%3N)
		#
	cd ..
	#
mkdir TEST_CSEM_MASSIF_INV
	#
	cd TEST_CSEM_MASSIF_INV
		#
		mpirun -np 2 valgrind --tool=massif --time-unit=ms --max-snapshots=1000 ../Mod3DMT -I NLCG ../../inputs/1st_Example/rFile_Model_trim ../../inputs/1st_Example/rFile_Data_CSEM ../../inputs/Naser_Zap/Inv_para.dat ../../inputs/Naser_Zap/FWD_para.dat | tee -a output.txt
		#
		bash ../../scripts/memory_brief.sh
		#
		T_END=$(date +%s%3N)
		echo "	> TEST_CSEM_MASSIF_INV: $(( ( $T_END - $T_START ) / 1000 )) seconds" | tee -a ../summary.txt
		#
	cd ..
	#
#
echo "MAIN SCRIPT FINISHES: $(date "+%Y/%m/%d - %H:%M:%S")" | tee -a summary.txt
#
exit 0
#
# END OF SCRIPT
