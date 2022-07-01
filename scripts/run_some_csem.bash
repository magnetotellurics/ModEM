#!/bin/bash
#
echo "MAIN SCRIPT STARTS: $(date "+%Y/%m/%d - %H:%M:%S")" | tee -a summary.txt
#
T_START=$(date +%s%3N)
#
# FWD CSEM
#
mkdir TEST_CSEM_FWD
	#
	cd TEST_CSEM_FWD
		#
		../test_FWD --forward --model ../../inputs/Naser_Zap/Start_model_WithGrad.mod --data ../../inputs/Naser_Zap/Field_CSEM_data20%err.dat --control ../../inputs/Others/first_control_file.txt | tee -a output.txt
		#
		T_END=$(date +%s%3N)
		echo "	> TEST_CSEM_FWD: $(( ( $T_END - $T_START ) / 1000 )) seconds" | tee -a ../summary.txt
		T_START=$(date +%s%3N)
		#
	cd ..
	#
mkdir TEST_CSEM_VALGRIND_FWD
	#
	cd TEST_CSEM_VALGRIND_FWD
		#
		valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --suppressions=../../../MEETING_INPUTS/valgring_mpi_suppression.supp --log-file=mt_valgrind_fwd_mpi.txt ../test_FWD --forward --model ../../inputs/Naser_Zap/Start_model_WithGrad.mod --data ../../inputs/Naser_Zap/Field_CSEM_data20%err.dat --control ../../inputs/Others/first_control_file.txt | tee -a output.txt
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
		valgrind --tool=massif --time-unit=ms --max-snapshots=1000 ../test_FWD --forward --model ../../inputs/Naser_Zap/Start_model_WithGrad.mod --data ../../inputs/Naser_Zap/Field_CSEM_data20%err.dat --control ../../inputs/Others/first_control_file.txt | tee -a output.txt
		#
		bash ../../scripts/memory_brief.sh
		#
		T_END=$(date +%s%3N)
		echo "	> TEST_CSEM_MASSIF_FWD: $(( ( $T_END - $T_START ) / 1000 )) seconds" | tee -a ../summary.txt
		T_START=$(date +%s%3N)
		#
	cd ..
	#
mkdir TEST_CSEM_FWD_MPI
	#
	cd TEST_CSEM_FWD_MPI
		#
		mpirun -np 7 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --suppressions=../../../MEETING_INPUTS/valgring_mpi_suppression.supp --log-file=mt_fwd_mpi.txt ../test_FWD_MPI --forward --model ../../inputs/Naser_Zap/Start_model_WithGrad.mod --data ../../inputs/Naser_Zap/Field_CSEM_data20%err.dat --control ../../inputs/Others/first_control_file.txt | tee -a output.txt
		#
		T_END=$(date +%s%3N)
		echo "	> TEST_CSEM_FWD_MPI: $(( ( $T_END - $T_START ) / 1000 )) seconds" | tee -a ../summary.txt
		T_START=$(date +%s%3N)
		#
	cd ..
	#
mkdir TEST_CSEM_VALGRIND_FWD_MPI
	#
	cd TEST_CSEM_VALGRIND_FWD_MPI
		#
		mpirun -np 7 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --suppressions=../../../MEETING_INPUTS/valgring_mpi_suppression.supp --log-file=mt_valgrind_fwd_mpi.txt ../test_FWD_MPI --forward --model ../../inputs/Naser_Zap/Start_model_WithGrad.mod --data ../../inputs/Naser_Zap/Field_CSEM_data20%err.dat --control ../../inputs/Others/first_control_file.txt | tee -a output.txt
		#
		T_END=$(date +%s%3N)
		echo "	> TEST_CSEM_VALGRIND_FWD_MPI: $(( ( $T_END - $T_START ) / 1000 )) seconds" | tee -a ../summary.txt
		T_START=$(date +%s%3N)
		#
	cd ..
	#
mkdir TEST_CSEM_MASSIF_FWD_MPI
	#
	cd TEST_CSEM_MASSIF_FWD_MPI
		#
		mpirun -np 7 valgrind --tool=massif --time-unit=ms --max-snapshots=1000 ../test_FWD_MPI --forward --model ../../inputs/Naser_Zap/Start_model_WithGrad.mod --data ../../inputs/Naser_Zap/Field_CSEM_data20%err.dat --control ../../inputs/Others/first_control_file.txt | tee -a output.txt
		#
		bash ../../scripts/memory_brief.sh
		#
		T_END=$(date +%s%3N)
		echo "	> TEST_CSEM_MASSIF_FWD_MPI: $(( ( $T_END - $T_START ) / 1000 )) seconds" | tee -a ../summary.txt
		T_START=$(date +%s%3N)
		#
	cd ..
	#
echo "MAIN SCRIPT FINISHES: $(date "+%Y/%m/%d - %H:%M:%S")" | tee -a summary.txt
#
exit 0
#
# END OF SCRIPT
