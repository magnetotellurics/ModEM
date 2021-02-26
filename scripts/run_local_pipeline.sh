#!/bin/#bash
#
#
NCORES=$1
#
#
NOW=$(date "+%Y_%m_%d_%H_%M_%S")
#
PIPELINE_TIME_START=$(date +%s%3N)
#
echo "OVER $(inputs/*)" | tee -a outputs/temp/summary.txt
#
# CREATE MAIN OUTPUT FOLDER (IF NOT EXISSTART)
mkdir -p outputs/temp
#
echo "$(hostname) START PIPELINE AT $NOW" | tee -a outputs/temp/summary.txt
echo "#" >> outputs/temp/summary.txt
echo "#" >> outputs/temp/summary.txt
#
# COMPILE MODEM STD AND SP2
echo "> START BUILD MODEM" >> outputs/temp/summary.txt
TIME_START=$(date +%s%3N)
bash scripts/build_modem_on.sh | tee outputs/temp/build_modem_on.txt
TIME_END=$(date +%s%3N)
echo "> FINISH BUILD MODEM" >> outputs/temp/summary.txt
echo "> Time Spent: $(( ( $TIME_END - $TIME_START ) / 1000 )) seconds" | tee -a outputs/temp/summary.txt
echo "#" >> outputs/temp/summary.txt
#
# TEST READ_WRITE STD
#
mkdir -p outputs/temp/test_read_write
#
# STD READ WRITE TEST WITH NPROC
echo "> START STD READ WRITE TEST" >> outputs/temp/summary.txt
TIME_START=$(date +%s%3N)
bash scripts/test_read_write.sh bin/ModEM_STD $NCORES | tee -a outputs/temp/test_read_write/std_out.txt
TIME_END=$(date +%s%3N)
echo "> FINISH STD READ WRITE TEST" >> outputs/temp/summary.txt
echo "> Time Spent: $(( ( $TIME_END - $TIME_START ) / 1000 )) seconds" | tee -a outputs/temp/summary.txt
echo "#" >> outputs/temp/summary.txt
#
# SP2 READ WRITE TEST WITH NPROC
echo "> START SP2 READ WRITE TEST" >> outputs/temp/summary.txt
TIME_START=$(date +%s%3N)
bash scripts/test_read_write.sh bin/ModEM_SP2 $NCORES | tee -a outputs/temp/test_read_write/std_out.txt
TIME_END=$(date +%s%3N)
echo "> FINISH SP2 READ WRITE TEST" >> outputs/temp/summary.txt
echo "> Time Spent: $(( ( $TIME_END - $TIME_START ) / 1000 )) seconds" | tee -a outputs/temp/summary.txt
echo "#" >> outputs/temp/summary.txt
#
# TEST FORWARD
mkdir -p outputs/temp/test_forward
#
# STD TEST FORWARD WITH NPROC
echo "> START STD FORWARD TEST" >> outputs/temp/summary.txt
TIME_START=$(date +%s%3N)
bash scripts/test_forward.sh bin/ModEM_STD $NCORES | tee -a outputs/temp/test_forward/std_out.txt
TIME_END=$(date +%s%3N)
echo "> FINISH STD FORWARD TEST" >> outputs/temp/summary.txt
echo "> Time Spent: $(( ( $TIME_END - $TIME_START ) / 1000 )) seconds" | tee -a outputs/temp/summary.txt
echo "#" >> outputs/temp/summary.txt
#
# SP2 TEST FORWARD WITH NPROC
echo "> START SP2 FORWARD TEST" >> outputs/temp/summary.txt
TIME_START=$(date +%s%3N)
bash scripts/test_forward.sh bin/ModEM_SP2 $NCORES | tee -a outputs/temp/test_forward/std_out.txt
TIME_END=$(date +%s%3N)
echo "> FINISH SP2 FORWARD TEST" >> outputs/temp/summary.txt
echo "> Time Spent: $(( ( $TIME_END - $TIME_START ) / 1000 )) seconds" | tee -a outputs/temp/summary.txt
echo "#" >> outputs/temp/summary.txt
#
# TEST MULT BY J T
mkdir -p outputs/temp/test_mult_by_j_t
#
# STD TEST MULT BY J T
echo "> START STD MULT BY J T TEST" >> outputs/temp/summary.txt
TIME_START=$(date +%s%3N)
bash scripts/test_mult_by_j_t.sh bin/ModEM_STD $NCORES | tee -a outputs/temp/test_mult_by_j_t/std_out.txt
#
RESULT_MULT_STD=$?
#
TIME_END=$(date +%s%3N)
echo "> FINISH STD MULT BY J T TEST" >> outputs/temp/summary.txt
echo "> Time Spent: $(( ( $TIME_END - $TIME_START ) / 1000 )) seconds" | tee -a outputs/temp/summary.txt
echo "#" >> outputs/temp/summary.txt
#
# SP2 TEST MULT BY J T
echo "> START SP2 MULT BY J T TEST" >> outputs/temp/summary.txt
TIME_START=$(date +%s%3N)
bash scripts/test_mult_by_j_t.sh bin/ModEM_SP2 $NCORES | tee -a outputs/temp/test_mult_by_j_t/std_out.txt
#
RESULT_MULT_SP2=$?
#
TIME_END=$(date +%s%3N)
echo "> FINISH SP2 MULT BY J T TEST" >> outputs/temp/summary.txt
echo "> Time Spent: $(( ( $TIME_END - $TIME_START ) / 1000 )) seconds" | tee -a outputs/temp/summary.txt
echo "#" >> outputs/temp/summary.txt
#
# TEST SYMMETRY
mkdir -p outputs/temp/test_symmetry
#
#if [ "$RESULT_MULT_STD" -ne "0" ]; then
	#
	# STD TEST SYMMETRY
	echo "> START STD SYMMETRY TEST" >> outputs/temp/summary.txt
	TIME_START=$(date +%s%3N)
	bash scripts/test_symmetry.sh bin/ModEM_STD $NCORES | tee -a outputs/temp/test_symmetry/std_out.txt
	TIME_END=$(date +%s%3N)
	echo "> FINISH STD SYMMETRY TEST" >> outputs/temp/summary.txt
	echo "> Time Spent: $(( ( $TIME_END - $TIME_START ) / 1000 )) seconds" | tee -a outputs/temp/summary.txt
	echo "#" >> outputs/temp/summary.txt
#else
	echo "> UNABLE TO START STD SYMMETRY TEST" >> outputs/temp/summary.txt
#fi
#
#
#if [ "$RESULT_MULT_SP2" -ne "0" ]; then
	#
	# SP2 TEST SYMMETRY
	echo "> START SP2 SYMMETRY TEST" >> outputs/temp/summary.txt
	TIME_START=$(date +%s%3N)
	bash scripts/test_symmetry.sh bin/ModEM_SP2 $NCORES | tee -a outputs/temp/test_symmetry/std_out.txt
	TIME_END=$(date +%s%3N)
	echo "> FINISH SP2 SYMMETRY TEST" >> outputs/temp/summary.txt
	echo "> Time Spent: $(( ( $TIME_END - $TIME_START ) / 1000 )) seconds" | tee -a outputs/temp/summary.txt
	echo "#" >> outputs/temp/summary.txt
#else
	echo "> UNABLE TO START SP2 SYMMETRY TEST" >> outputs/temp/summary.txt
#fi
#
# HANDLE OUTPUT
rm -rf bin/
#
PIPELINE_TIME_END=$(date +%s%3N)
echo "#" >> outputs/temp/summary.txt
echo "$(hostname) FINISH PIPELINE" | tee -a outputs/temp/summary.txt
echo "> Total Time Spent: $(( ( $PIPELINE_TIME_END - $PIPELINE_TIME_START ) / 1000 )) seconds" | tee -a outputs/temp/summary.txt
#
# CREATE FOLDER UNIQUE NAME
mkdir -p outputs/$(hostname)_pipeline_$NOW
#
# CHANGE TEMP NAME
cp -r outputs/temp/* outputs/$(hostname)_pipeline_$NOW/
#
rm -rf outputs/temp
#
# END MAIN SCRIPT