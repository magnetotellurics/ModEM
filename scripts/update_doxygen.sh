#!/bin/bash
#
# 
is_local=$1
#
echo "#### Start Doxygen documentation at $(date "+%Y/%m/%d - %H:%M:%S")" | tee -a outputs/temp/summary.txt
echo "    #" | tee -a outputs/temp/summary.txt
#
T_START=$(date +%s%3N)
#
cd src
    #
    # CLEAN
    doxygen ../docs/doxygen_modem_oo_config
    #
    # CATCH RESULT
    RESULT=$?
    #
    T_END=$(date +%s%3N)
    #
    # TEST RESULT
    if [ "$RESULT" -ne "0" ]; then
        #
        echo "    > Doxygen documentation fail: $RESULT" | tee -a ../outputs/temp/summary.txt
        echo "    > Time Spent: $(( ( $T_END - $T_START ) / 1000 )) seconds" | tee -a ../outputs/temp/summary.txt
        echo "    #" | tee -a ../outputs/temp/summary.txt
        #
        cd ..
        #
        exit $RESULT
        #
    fi
    #
    echo "    > Doxygen documentation successful" | tee -a ../outputs/temp/summary.txt
    echo "    > Time Spent: $(( ( $T_END - $T_START ) / 1000 )) seconds" | tee -a ../outputs/temp/summary.txt
    echo "    #" | tee -a ../outputs/temp/summary.txt
    #
    cd ..
#
echo "#### Finish build ModEM-OO at $(date "+%Y/%m/%d - %H:%M:%S")" | tee -a outputs/temp/summary.txt
#
if $is_local ; then
	echo 'is local'
	# RENAME outputs/temp
	mv outputs/temp outputs/$(hostname)_$(date "+%Y%m%d_%H%M%S")
	#
fi
#
exit 0
#
# END OF SCRIPT

