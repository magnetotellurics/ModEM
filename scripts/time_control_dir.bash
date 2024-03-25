#!/bin/bash
#
target_dir=$1
#
procs=9
#
FIRST=1
LAST=10
#
start_model="MR_TEST_INPUTS/MT/m0.rho"
real_model="MR_TEST_INPUTS/MT/m.rho"
data="MR_TEST_INPUTS/MT/d_05.dat"
#
FWD_CONTROLS=$(ls ${target_dir}*)
#
for FWD_CONTROL in $FWD_CONTROLS; do
    #
    #echo $FWD_CONTROL
    #
    OUTPUT=$FWD_CONTROL
    OUTPUT=${OUTPUT/.exe/$emptyspace}
    OUTPUT=${OUTPUT/.sh/$emptyspace}
    OUTPUT=${OUTPUT/.txt/$emptyspace}
    OUTPUT=${OUTPUT/*\//$emptyspace}
    #
    T_START=$(date +%s%3N)
    #
    for (( i=$FIRST; i<=$LAST; i++ )); do
        #
        #echo $i
        #
        mpirun -np $procs ./ModEM -f -m $real_model -d $data -cf $FWD_CONTROL -pd PRED_$OUTPUT.dat | tee -a PRED_$OUTPUT.txt
        #
        #mpirun -np 9 ./Mod3DMT_STD -F ../../modem-oo/src/MR_TEST_INPUTS/MT/m.rho ../../modem-oo/src/MR_TEST_INPUTS/MT/d_05.dat PRED_ON.dat esol.bin FWD_para.dat
        #
    done
    #
    T_END=$(date +%s%3N)
    #
    echo $OUTPUT $(( ( $T_END - $T_START ) / ( 1000 * $LAST ) )) | tee -a FWD_TIMES.txt
    #
done
#
exit 0
#
# LAST OF SCRIPT
