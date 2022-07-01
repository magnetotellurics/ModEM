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
    mkdir TEST_MT_FWD_MPI
        #
        cd TEST_MT_FWD_MPI
            #
            mpirun -np 2 ../../Mod3DMT -F ../../../../MEETING_INPUTS/rFile_Model ../../../../MEETING_INPUTS/rFile_Data_MT wFile_Data wFile_EMsoln ../../../../MEETING_INPUTS/rFile_fwdCtrl | tee -a output.txt
            #
            T_END=$(date +%s%3N)
            echo "    > TEST_MT_FWD_MPI: $(( ( $T_END - $T_START ) / 1000 )) seconds" | tee -a ../summary.txt
            T_START=$(date +%s%3N)
            #
        cd ..
        #
    mkdir TEST_MT_VALGRIND_FWD_MPI
        #
        cd TEST_MT_VALGRIND_FWD_MPI
            #
            mpirun -np 2 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --suppressions=../../../../MEETING_INPUTS/valgring_mpi_suppression.supp --log-file=mt_valgrind_fwd_mpi.txt ../../Mod3DMT -F ../../../../MEETING_INPUTS/rFile_Model ../../../../MEETING_INPUTS/rFile_Data_MT wFile_Data wFile_EMsoln ../../../../MEETING_INPUTS/rFile_fwdCtrl | tee -a output.txt
            #
            T_END=$(date +%s%3N)
            echo "    > TEST_MT_VALGRIND_FWD_MPI: $(( ( $T_END - $T_START ) / 1000 )) seconds" | tee -a ../summary.txt
            T_START=$(date +%s%3N)
            #
        cd ..
        #
    mkdir TEST_MT_MASSIF_FWD_MPI
        #
        cd TEST_MT_MASSIF_FWD_MPI
            #
            mpirun -np 2 valgrind --tool=massif --time-unit=ms --max-snapshots=1000 ../../Mod3DMT -F ../../../../MEETING_INPUTS/rFile_Model ../../../../MEETING_INPUTS/rFile_Data_MT wFile_Data wFile_EMsoln ../../../../MEETING_INPUTS/rFile_fwdCtrl | tee -a output.txt
            #
            bash ../../scripts/memory_brief.sh
            #
            T_END=$(date +%s%3N)
            echo "    > TEST_MT_MASSIF_FWD_MPI: $(( ( $T_END - $T_START ) / 1000 )) seconds" | tee -a ../summary.txt
            T_START=$(date +%s%3N)
            #
        cd ..
    #
    # FWD CSEM
    #
    mkdir TEST_CSEM_FWD_MPI
        #
        cd TEST_CSEM_FWD_MPI
            #
            mpirun -np 2 ../../Mod3DMT -F ../../../../MEETING_INPUTS/rFile_Model ../../../../MEETING_INPUTS/rFile_Data_CSEM wFile_Data wFile_EMsoln ../../../../MEETING_INPUTS/rFile_fwdCtrl | tee -a output.txt
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
            mpirun -np 2 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --suppressions=../../../../MEETING_INPUTS/valgring_mpi_suppression.supp --log-file=mt_valgrind_fwd_mpi.txt ../../Mod3DMT -F ../../../../MEETING_INPUTS/rFile_Model ../../../../MEETING_INPUTS/rFile_Data_CSEM wFile_Data wFile_EMsoln ../../../../MEETING_INPUTS/rFile_fwdCtrl | tee -a output.txt
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
            mpirun -np 2 valgrind --tool=massif --time-unit=ms --max-snapshots=1000 ../../Mod3DMT -F ../../../../MEETING_INPUTS/rFile_Model ../../../../MEETING_INPUTS/rFile_Data_CSEM wFile_Data wFile_EMsoln ../../../../MEETING_INPUTS/rFile_fwdCtrl | tee -a output.txt
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
    # INV MT
    #
    mkdir TEST_MT_INV_MPI
        #
        cd TEST_MT_INV_MPI
            #
            mpirun -np 2 ../../Mod3DMT -I NLCG ../../../../MEETING_INPUTS/rFile_Model ../../../../MEETING_INPUTS/rFile_Data_MT ../../../../MEETING_INPUTS/rFile_invCtrl ../../../../MEETING_INPUTS/rFile_fwdCtrl | tee -a output.txt
            #
            T_END=$(date +%s%3N)
            echo "    > TEST_MT_INV_MPI: $(( ( $T_END - $T_START ) / 1000 )) seconds" | tee -a ../summary.txt
            T_START=$(date +%s%3N)
            #
        cd ..
        #
    mkdir TEST_MT_VALGRIND_INV_MPI
        #
        cd TEST_MT_VALGRIND_INV_MPI
            #
            mpirun -np 2 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --suppressions=../../../../MEETING_INPUTS/valgring_mpi_suppression.supp --log-file=mt_valgrind_fwd_mpi.txt ../../Mod3DMT -I NLCG ../../../../MEETING_INPUTS/rFile_Model ../../../../MEETING_INPUTS/rFile_Data_MT ../../../../MEETING_INPUTS/rFile_invCtrl ../../../../MEETING_INPUTS/rFile_fwdCtrl | tee -a output.txt
            #
            T_END=$(date +%s%3N)
            echo "    > TEST_MT_VALGRIND_INV_MPI: $(( ( $T_END - $T_START ) / 1000 )) seconds" | tee -a ../summary.txt
            T_START=$(date +%s%3N)
            #
        cd ..
        #
    mkdir TEST_MT_MASSIF_INV_MPI
        #
        cd TEST_MT_MASSIF_INV_MPI
            #
            mpirun -np 2 valgrind --tool=massif --time-unit=ms --max-snapshots=1000 ../../Mod3DMT -I NLCG ../../../../MEETING_INPUTS/rFile_Model ../../../../MEETING_INPUTS/rFile_Data_MT ../../../../MEETING_INPUTS/rFile_invCtrl ../../../../MEETING_INPUTS/rFile_fwdCtrl | tee -a output.txt
            #
            bash ../../scripts/memory_brief.sh
            #
            T_END=$(date +%s%3N)
            echo "    > TEST_MT_MASSIF_INV_MPI: $(( ( $T_END - $T_START ) / 1000 )) seconds" | tee -a ../summary.txt
            T_START=$(date +%s%3N)
            #
        cd ..
    #
    # INV CSEM
    #
    mkdir TEST_CSEM_INV_MPI
        #
        cd TEST_CSEM_INV_MPI
            #
            mpirun -np 2 ../../Mod3DMT -I NLCG ../../../../MEETING_INPUTS/rFile_Model ../../../../MEETING_INPUTS/rFile_Data_CSEM ../../../../MEETING_INPUTS/rFile_invCtrl ../../../../MEETING_INPUTS/rFile_fwdCtrl | tee -a output.txt
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
            mpirun -np 2 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --suppressions=../../../../MEETING_INPUTS/valgring_mpi_suppression.supp --log-file=mt_valgrind_fwd_mpi.txt ../../Mod3DMT -I NLCG ../../../../MEETING_INPUTS/rFile_Model ../../../../MEETING_INPUTS/rFile_Data_CSEM ../../../../MEETING_INPUTS/rFile_invCtrl ../../../../MEETING_INPUTS/rFile_fwdCtrl | tee -a output.txt
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
            mpirun -np 2 valgrind --tool=massif --time-unit=ms --max-snapshots=1000 ../../Mod3DMT -I NLCG ../../../../MEETING_INPUTS/rFile_Model ../../../../MEETING_INPUTS/rFile_Data_CSEM ../../../../MEETING_INPUTS/rFile_invCtrl ../../../../MEETING_INPUTS/rFile_fwdCtrl | tee -a output.txt
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
