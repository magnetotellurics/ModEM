#!/bin/bash
#
procs=9
start_model="MR_TEST_INPUTS/MT/m0.rho"
real_model="MR_TEST_INPUTS/MT/m.rho"
data="MR_TEST_INPUTS/MT/d05.dat"
#
FWD_CONTROLS=$(ls CONTROL_FILES/*)
#
for FWD_CONTROL in $FWD_CONTROLS; do
    #
    echo $FWD_CONTROL
    #
    OUTPUT=$FWD_CONTROL
    OUTPUT=${OUTPUT/.exe/$emptyspace}
    OUTPUT=${OUTPUT/.sh/$emptyspace}
    OUTPUT=${OUTPUT/.txt/$emptyspace}
    OUTPUT=${OUTPUT/*\//$emptyspace}
    #
    T_START=$(date +%s%3N)
    # #
    # valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=fwd_$OUTPUT.vg ./ModEM_SERIAL -f -m $real_model -d $data -cf $FWD_CONTROL -pd PRED_$OUTPUT.dat | tee -a PRED_$OUTPUT.txt
    # #
    # valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=jm_$OUTPUT.vg ./ModEM_SERIAL -j -m $start_model -pm $real_model -d $data -cf $FWD_CONTROL -jm JHMAT_$OUTPUT.dat | tee -a JHMAT_$OUTPUT.txt
    # #
    # valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=jt_$OUTPUT.vg ./ModEM_SERIAL -jt -m $start_model -d $data -cf $FWD_CONTROL -dm DSIGMA_$OUTPUT.rho | tee -a DSIGMA_$OUTPUT.txt
    # #
    mpirun -np $procs ./ModEM_MPI -f -m $real_model -d $data -cf $FWD_CONTROL -pd PRED_$OUTPUT.dat | tee -a PRED_$OUTPUT.txt
    #
    mpirun -np $procs ./ModEM_MPI -j -m $start_model -pm $real_model -d $data -cf $FWD_CONTROL -jm JHMAT_$OUTPUT.dat | tee -a JHMAT_$OUTPUT.txt
    #
    mpirun -np $procs ./ModEM_MPI -jt -m $start_model -d $data -cf $FWD_CONTROL -dm DSIGMA_$OUTPUT.rho | tee -a DSIGMA_$OUTPUT.txt
    #
    T_END=$(date +%s%3N)
    #
    echo $OUTPUT $(( ( $T_END - $T_START ) / 1000 ))
    #
done
#
exit 0
#
# END OF SCRIPT