#!/bin/bash
#
procs=9
#
FWD_CONTROLS=$(ls CONTROL_FILES/*)
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
    #for i in {1..100}; do 
        #
        #echo $i
        mpirun -np $procs valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v --log-file=$OUTPUT.valg ./ModEM_MPI -f -m MR_TEST_INPUTS/m.rho -d MR_TEST_INPUTS/d0.dat -cf $FWD_CONTROL -pd $OUTPUT.dat > $OUTPUT.txt
        #
    #done
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

