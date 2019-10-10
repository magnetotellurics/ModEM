#!/bin/bash
#
# ARGUMENTS: 1 - MODEM PATH
#
EXEC=$1
#
#
emptyspace=""
#
# STRING NOW
now=$(date "+%Y/%m/%d - %H:%M:%S")
#
#
# FOREACH TEST DIR ON tests/
for TEST_DIR in tests/*
do
	cd $TEST_DIR
	#
	# VAR DIR_NAME
	test_name=$TEST_DIR
	test_name=${test_name/*\//$emptyspace}
	#
	# EXECUTE BASH RUN 
	bash run.sh ../$EXEC | tee std_out.txt
	#
	# CATCH RESULT
    result=$?
	#
	# TEST RESULT
	if [ "$result" -ne "0" ]; then
	   echo "$test_name FAIL: $result" | tee -a std_out.txt
	   exit $result
	fi
	#
	echo "$test_name PASS" | tee -a std_out.txt
	#
	# REMOVE OLD TEST FOLDER RESULT FROM outputs/
	rm -rf ../../outputs/result_${test_name}
	#
	# CREATE NEW TEST FOLDER RESULT ON outputs/
	mkdir ../../outputs/result_${test_name}
	#
	# MOVE ALL FILES TO NEW TEST FOLDER RESULT ON outputs/
	mv * ../../outputs/result_${test_name}
	#
	#
	cd ..
done
#
exit 0
#
#
# END OF SCRIPT
