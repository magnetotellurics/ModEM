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
now=$(date "+%Y_%m_%d_%H_%M_%S")
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
	# CREATE OUTPUT FOLDER FOR THE TEST
	mkdir -p results_${test_name}_$now
	#
	# EXECUTE BASH RUN 
	bash run.sh ../$EXEC &>> ${test_name}_std_out.txt
	#
	# CATCH RESULT
    	result=$?
	#
	# TEST RESULT
	if [ "$result" -ne "0" ]; then
	   echo "$test_name FAIL: $result" &>> ${test_name}_std_out.txt
	   exit $result
	fi
	#
	echo "$test_name PASS" &>> ${test_name}_std_out.txt
	#
	# MOVES EVERYTHING STARTING WITH THE TEST NAME TO OUTPUT FOLDER
	mv ${test_name}* results_${test_name}_$now/
	#
	# MOVES TEST OUTPUT FOLDER TO MAIN OUTPUT FOLDER
	mv results_${test_name}_$now/ ../../outputs/
	#
	cd ..
	#
	cd ..
done
#
exit 0
#
#
# END OF SCRIPT
