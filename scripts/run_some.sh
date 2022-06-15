#!/bin/bash
#
echo "MAIN SCRIPT STARTS: $(date "+%Y/%m/%d - %H:%M:%S")" | tee -a summary.txt
#
T_START=$(date +%s%3N)
#
# FWD
#
mkdir TEST_MT_VALGRIND_FWD
	#
	cd TEST_MT_VALGRIND_FWD
		#
		mpirun -np 9  valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --log-file=fwd_valgrind_mt.txt ../Mod3DMT -F ../../inputs/Benchmarking/rFile_Model ../../inputs/Benchmarking/rFile_Data predicted
		#
		T_END=$(date +%s%3N)
		echo "	> TEST_MT_VALGRIND_FWD: $(( ( $T_END - $T_START ) / 1000 )) seconds" | tee -a ../summary.txt
		T_START=$(date +%s%3N)
		#
	cd ..
	#
mkdir TEST_CSEM_VALGRIND_FWD
	#
	cd TEST_CSEM_VALGRIND_FWD
		#
		mpirun -np 46  valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --log-file=fwd_valgrind_csem.txt ../Mod3DMT -F ../../inputs/Naser_Zap/Start_model_WithGrad.mod  ../../inputs/Naser_Zap/Field_CSEM_data20%err.dat predicted
		#
		T_END=$(date +%s%3N)
		echo "	> TEST_CSEM_VALGRIND_FWD: $(( ( $T_END - $T_START ) / 1000 )) seconds" | tee -a ../summary.txt
		T_START=$(date +%s%3N)
		#
	cd ..
	#
	#
mkdir TEST_MT_MASSIF_FWD
	#
	cd TEST_MT_MASSIF_FWD
		#
		mpirun -np 9 valgrind --tool=massif --time-unit=ms --max-snapshots=1000 ../Mod3DMT -F ../../inputs/Benchmarking/rFile_Model ../../inputs/Benchmarking/rFile_Data predicted
		#
		bash ../../scripts/memory_brief.sh
		#
		T_END=$(date +%s%3N)
		echo "	> TEST_MT_MASSIF_FWD: $(( ( $T_END - $T_START ) / 1000 )) seconds" | tee -a ../summary.txt
		T_START=$(date +%s%3N)
		#
	cd ..
	#
mkdir TEST_CSEM_MASSIF_FWD
	#
	cd TEST_CSEM_MASSIF_FWD
		#
		mpirun -np 46 valgrind --tool=massif --time-unit=ms --max-snapshots=1000 ../Mod3DMT -F ../../inputs/Naser_Zap/Start_model_WithGrad.mod  ../../inputs/Naser_Zap/Field_CSEM_data20%err.dat predicted
		#
		bash ../../scripts/memory_brief.sh
		#
		T_END=$(date +%s%3N)
		echo "	> TEST_CSEM_MASSIF_FWD: $(( ( $T_END - $T_START ) / 1000 )) seconds" | tee -a ../summary.txt
		T_START=$(date +%s%3N)
		#
	cd ..
	#
	#
#
# INV
#
mkdir TEST_MT_VALGRIND_INV
	#
	cd TEST_MT_VALGRIND_INV
		#
		mpirun -np 9  valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --log-file=inv_valgrind_mt.txt ../Mod3DMT -I NLCG ../../inputs/Benchmarking/rFile_Model ../../inputs/Benchmarking/rFile_Data
		#
		T_END=$(date +%s%3N)
		echo "	> TEST_MT_VALGRIND_INV: $(( ( $T_END - $T_START ) / 1000 )) seconds" | tee -a ../summary.txt
		T_START=$(date +%s%3N)
		#
	cd ..
	#
	#
mkdir TEST_CSEM_VALGRIND_INV
	#
	cd TEST_CSEM_VALGRIND_INV
		#
		mpirun -np 46  valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --log-file=inv_valgrind_csem.txt ../Mod3DMT -I NLCG ../../inputs/Naser_Zap/Start_model_WithGrad.mod ../../inputs/Naser_Zap/Field_CSEM_data20%err.dat ../../inputs/Naser_Zap/Inv_para.dat ../../inputs/Naser_Zap/FWD_para.dat ../../inputs/Naser_Zap/Start_model_WithGrad.cov
		#
		T_END=$(date +%s%3N)
		echo "	> TEST_CSEM_VALGRIND_INV: $(( ( $T_END - $T_START ) / 1000 )) seconds" | tee -a ../summary.txt
		T_START=$(date +%s%3N)
		#
	cd ..
	#
	#
mkdir TEST_MT_MASSIF_INV
	#
	cd TEST_MT_MASSIF_INV
		#
		mpirun -np 9 valgrind --tool=massif --time-unit=ms --max-snapshots=1000 ../Mod3DMT -I NLCG ../../inputs/Benchmarking/rFile_Model ../../inputs/Benchmarking/rFile_Data
		#
		bash ../../scripts/memory_brief.sh
		#
		T_END=$(date +%s%3N)
		echo "	> TEST_MT_MASSIF_INV: $(( ( $T_END - $T_START ) / 1000 )) seconds" | tee -a ../summary.txt
		T_START=$(date +%s%3N)
		#
	cd ..
	#
mkdir TEST_CSEM_MASSIF_INV
	#
	cd TEST_CSEM_MASSIF_INV
	#
	mpirun -np 46 valgrind --tool=massif --time-unit=ms --max-snapshots=1000 ../Mod3DMT -I NLCG ../../inputs/Naser_Zap/Start_model_WithGrad.mod ../../inputs/Naser_Zap/Field_CSEM_data20%err.dat ../../inputs/Naser_Zap/Inv_para.dat ../../inputs/Naser_Zap/FWD_para.dat ../../inputs/Naser_Zap/Start_model_WithGrad.cov
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
