#!/bin/bash
#
target_dir=$1
#
procs=9
#
FIRST=1
LAST=1
#
start_model="../inputs/MR_TEST_INPUTS/MT/m0.rho"
real_model="../inputs/MR_TEST_INPUTS/MT/m.rho"
data="../inputs/MR_TEST_INPUTS/MT/d_05.dat"
#
FWD_CONTROLS=$(ls ${target_dir}*)
#
for FWD_CONTROL in $FWD_CONTROLS; do
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
        #mpirun -np $procs ./ModEM_MPI -f -m $real_model -d $data -cf $FWD_CONTROL -pd PRED_$OUTPUT.dat | tee -a PRED_$OUTPUT.txt
        #
        #mpirun -np $procs ./ModEM_MPI -j -m $start_model -pm $real_model -d $data -cf $FWD_CONTROL -jm JHMAT_$OUTPUT.dat | tee -a JHMAT_$OUTPUT.txt
        #
        mpirun -np $procs ./ModEM_MPI -jt -m $start_model -d $data -cf $FWD_CONTROL -dm DSIGMA_$OUTPUT.rho | tee -a DSIGMA_$OUTPUT.txt
        #
        #mpirun -np $procs ./ModEM_MPI -i -m $start_model -d $data -ci inv_ctrl_template.txt -cf $FWD_CONTROL -o INV_$OUTPUT.rho | tee -a INV_$OUTPUT.txt
        #
    done
    #
    T_END=$(date +%s%3N)
    #
    echo $OUTPUT $(( ( $T_END - $T_START ) / ( 1000 * $LAST ) )) | tee -a JT_TIMES.txt
    #
done
#
exit 0
#
# END OF SCRIPT