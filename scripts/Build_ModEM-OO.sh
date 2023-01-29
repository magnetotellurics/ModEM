#!/bin/bash
#
is_local=$1
#
# CREATE BIN FOLDER (IF NOT EXISTS)
mkdir -p bin
#
# CREATE MAIN OUTPUT FOLDER (IF NOT EXISTS)
mkdir -p outputs/temp
#
echo "#### Start build ModEM-OO at $(date "+%Y/%m/%d - %H:%M:%S")" | tee -a outputs/temp/summary.txt
echo "    #" | tee -a outputs/temp/summary.txt
#
T_START=$(date +%s%3N)
#
# REMOVE OLD OBJECTS
OBJ_DIR="../objs"
#
if [ -d "$DIR" ]; then
    rm -rf $OBJ_DIR
fi
#
#
cd src
    #
    # CLEAN
    make -f MakefileSerial clean | tee -a ../outputs/temp/build_output.txt
    #
    # BUILD MAKEFILE
    ./Configure.modem.oo.GFortran MakefileSerial SERIAL TestSerial.f90 | tee -a ../outputs/temp/build_output.txt
    #
    # BUILD EXECUTABLE WITH THE MAKEFILE
    make -f MakefileSerial | tee -a ../outputs/temp/build_output.txt
    #
    # CATCH RESULT
    RESULT=$?
    #
    T_END=$(date +%s%3N)
    #
    # TEST RESULT
    if [ "$RESULT" -ne "0" ]; then
        #
        echo "    > Build ModEM-OO fail: $RESULT" | tee -a ../outputs/temp/summary.txt
        echo "    > Time Spent: $(( ( $T_END - $T_START ) / 1000 )) seconds" | tee -a ../outputs/temp/summary.txt
        echo "    #" | tee -a ../outputs/temp/summary.txt
        #
        cd ..
        #
        exit $RESULT
        #
    fi
    #
    echo "    > Build ModEM-OO successful" | tee -a ../outputs/temp/summary.txt
    echo "    > Time Spent: $(( ( $T_END - $T_START ) / 1000 )) seconds" | tee -a ../outputs/temp/summary.txt
    echo "    #" | tee -a ../outputs/temp/summary.txt
    #
    # RENAME EXE
    mv TestSerial ../bin
    #
    cd ..
#
# REMOVE NEW OBJECTS
rm -rf $OBJ_DIR
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

