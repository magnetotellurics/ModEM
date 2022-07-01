#!/bin/bash
#
OUTPUT_FOLDER=$1
#
T_START=$(date +%s%3N)
#
# FWD MT
#
mkdir $OUTPUT_FOLDER
cd $OUTPUT_FOLDER
    #
    echo "MAIN SCRIPT STARTS: $(date "+%Y/%m/%d - %H:%M:%S")" | tee -a summary.txt
    #
    # FWD CSEM
    #
    mkdir TEST_CSEM_FWD_MPI
        #
        cd TEST_CSEM_FWD_MPI
            #
            mpirun -np 46 ../../Mod3DMT -F ../../../inputs/Naser_Zap/Start_model_WithGrad.mod ../../../inputs/Naser_Zap/Field_CSEM_data20%err.dat wFile_Data wFile_EMsoln ../../../../MEETING_INPUTS/rFile_fwdCtrl | tee -a output.txt
            #
            T_END=$(date +%s%3N)
            echo "    > TEST_CSEM_FWD_MPI: $(( ( $T_END - $T_START ) / 1000 )) seconds" | tee -a ../summary.txt
            T_START=$(date +%s%3N)
            #
        cd ..
        #
    mkdir TEST_CSEM_VALGRIND_FWD_MPI
        #
        cd TEST_CSEM_VALGRIND_FWD_MPI
            #
            mpirun -np 46 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --suppressions=../../../../MEETING_INPUTS/valgring_mpi_suppression.supp --log-file=mt_valgrind_fwd_mpi.txt ../../Mod3DMT -F ../../../inputs/Naser_Zap/Start_model_WithGrad.mod ../../../inputs/Naser_Zap/Field_CSEM_data20%err.dat wFile_Data wFile_EMsoln ../../../../MEETING_INPUTS/rFile_fwdCtrl | tee -a output.txt
            #
            T_END=$(date +%s%3N)
            echo "    > TEST_CSEM_VALGRIND_FWD_MPI: $(( ( $T_END - $T_START ) / 1000 )) seconds" | tee -a ../summary.txt
            T_START=$(date +%s%3N)
            #
        cd ..
        #
    mkdir TEST_CSEM_MASSIF_FWD_MPI
        #
        cd TEST_CSEM_MASSIF_FWD_MPI
            #
            mpirun -np 46 valgrind --tool=massif --time-unit=ms --max-snapshots=1000 ../../Mod3DMT -F ../../../inputs/Naser_Zap/Start_model_WithGrad.mod ../../../inputs/Naser_Zap/Field_CSEM_data20%err.dat wFile_Data wFile_EMsoln ../../../../MEETING_INPUTS/rFile_fwdCtrl | tee -a output.txt
            #
            bash ../../scripts/memory_brief.sh
            #
            T_END=$(date +%s%3N)
            echo "    > TEST_CSEM_MASSIF_FWD_MPI: $(( ( $T_END - $T_START ) / 1000 )) seconds" | tee -a ../summary.txt
            T_START=$(date +%s%3N)
            #
        cd ..
        #
    echo "FWD FINISHES: $(date "+%Y/%m/%d - %H:%M:%S")" | tee -a summary.txt
    #
    # INV CSEM
    #
    mkdir TEST_CSEM_INV_MPI
        #
        cd TEST_CSEM_INV_MPI
            #
            mpirun -np 46 ../../Mod3DMT -I NLCG ../../../inputs/Naser_Zap/Start_model_WithGrad.mod ../../../inputs/Naser_Zap/Field_CSEM_data20%err.dat ../../../../MEETING_INPUTS/rFile_invCtrl ../../../../MEETING_INPUTS/rFile_fwdCtrl | tee -a output.txt
            #
            T_END=$(date +%s%3N)
            echo "    > TEST_CSEM_INV_MPI: $(( ( $T_END - $T_START ) / 1000 )) seconds" | tee -a ../summary.txt
            T_START=$(date +%s%3N)
            #
        cd ..
        #
    mkdir TEST_CSEM_VALGRIND_INV_MPI
        #
        cd TEST_CSEM_VALGRIND_INV_MPI
            #
            mpirun -np 46 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --suppressions=../../../../MEETING_INPUTS/valgring_mpi_suppression.supp --log-file=mt_valgrind_fwd_mpi.txt ../../Mod3DMT -I NLCG ../../../inputs/Naser_Zap/Start_model_WithGrad.mod ../../../inputs/Naser_Zap/Field_CSEM_data20%err.dat ../../../../MEETING_INPUTS/rFile_invCtrl ../../../../MEETING_INPUTS/rFile_fwdCtrl | tee -a output.txt
            #
            T_END=$(date +%s%3N)
            echo "    > TEST_CSEM_VALGRIND_INV_MPI: $(( ( $T_END - $T_START ) / 1000 )) seconds" | tee -a ../summary.txt
            T_START=$(date +%s%3N)
            #
        cd ..
        #
    mkdir TEST_CSEM_MASSIF_INV_MPI
        #
        cd TEST_CSEM_MASSIF_INV_MPI
            #
            mpirun -np 46 valgrind --tool=massif --time-unit=ms --max-snapshots=1000 ../../Mod3DMT -I NLCG ../../../inputs/Naser_Zap/Start_model_WithGrad.mod ../../../inputs/Naser_Zap/Field_CSEM_data20%err.dat ../../../../MEETING_INPUTS/rFile_invCtrl ../../../../MEETING_INPUTS/rFile_fwdCtrl | tee -a output.txt
            #
            bash ../../scripts/memory_brief.sh
            #
            T_END=$(date +%s%3N)
            echo "    > TEST_CSEM_MASSIF_INV_MPI: $(( ( $T_END - $T_START ) / 1000 )) seconds" | tee -a ../summary.txt
            T_START=$(date +%s%3N)
            #
        cd ..
        #
    echo "MAIN SCRIPT FINISHES: $(date "+%Y/%m/%d - %H:%M:%S")" | tee -a summary.txt
    #
cd ..
exit 0
#
# END OF SCRIPT
