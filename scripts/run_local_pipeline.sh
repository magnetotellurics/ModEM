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
# CREATE MAIN OUTPUT FOLDER (IF NOT EXIST)
mkdir -p outputs/temp
#
echo "$(hostname) START PIPELINE AT $NOW" | tee -a outputs/temp/summary.txt
echo "#" >> outputs/temp/summary.txt
echo "#" >> outputs/temp/summary.txt
#
# COMPILE ModEM_Baseline
echo "> START BUILD MODEM" >> outputs/temp/summary.txt
TIME_START=$(date +%s%3N)
bash scripts/build_modem_on.sh | tee outputs/temp/build_modem_on.txt
TIME_END=$(date +%s%3N)
echo "> FINISH BUILD MODEM" >> outputs/temp/summary.txt
echo "> Time Spent: $(( ( $TIME_END - $TIME_START ) / 1000 )) seconds" | tee -a outputs/temp/summary.txt
echo "#" >> outputs/temp/summary.txt
#
# TEST READ_WRITE BASELINE
#
mkdir -p outputs/temp/test_read_write
#
# BASELINE READ WRITE TEST WITH NPROC
echo "> START BASELINE READ WRITE TEST" >> outputs/temp/summary.txt
TIME_START=$(date +%s%3N)
bash scripts/auto/test_read_write.sh bin/ModEM_Baseline $NCORES | tee -a outputs/temp/test_read_write/std_out.txt
TIME_END=$(date +%s%3N)
echo "> FINISH BASELINE READ WRITE TEST" >> outputs/temp/summary.txt
echo "> Time Spent: $(( ( $TIME_END - $TIME_START ) / 1000 )) seconds" | tee -a outputs/temp/summary.txt
echo "#" >> outputs/temp/summary.txt
#
# TEST FORWARD
mkdir -p outputs/temp/test_forward
#
# BASELINE TEST FORWARD WITH NPROC
echo "> START BASELINE FORWARD TEST" >> outputs/temp/summary.txt
TIME_START=$(date +%s%3N)
bash scripts/auto/test_forward.sh bin/ModEM_Baseline $NCORES | tee -a outputs/temp/test_forward/std_out.txt
TIME_END=$(date +%s%3N)
echo "> FINISH BASELINE FORWARD TEST" >> outputs/temp/summary.txt
echo "> Time Spent: $(( ( $TIME_END - $TIME_START ) / 1000 )) seconds" | tee -a outputs/temp/summary.txt
echo "#" >> outputs/temp/summary.txt
#
# TEST MULT BY J T
mkdir -p outputs/temp/test_mult_by_j_t
#
# BASELINE TEST MULT BY J T
echo "> START BASELINE MULT BY J T TEST" >> outputs/temp/summary.txt
TIME_START=$(date +%s%3N)
bash scripts/auto/test_mult_by_j_t.sh bin/ModEM_Baseline $NCORES | tee -a outputs/temp/test_mult_by_j_t/std_out.txt
#
RESULT_MULT_BASELINE=$?
#
TIME_END=$(date +%s%3N)
echo "> FINISH BASELINE MULT BY J T TEST" >> outputs/temp/summary.txt
echo "> Time Spent: $(( ( $TIME_END - $TIME_START ) / 1000 )) seconds" | tee -a outputs/temp/summary.txt
echo "#" >> outputs/temp/summary.txt
#
# TEST SYMMETRY
mkdir -p outputs/temp/test_symmetry
#
#if [ "$RESULT_MULT_BASELINE" -ne "0" ]; then
	#
	# BASELINE TEST SYMMETRY
	echo "> START BASELINE SYMMETRY TEST" >> outputs/temp/summary.txt
	TIME_START=$(date +%s%3N)
	bash scripts/auto/test_symmetry.sh bin/ModEM_Baseline $NCORES | tee -a outputs/temp/test_symmetry/std_out.txt
	TIME_END=$(date +%s%3N)
	echo "> FINISH BASELINE SYMMETRY TEST" >> outputs/temp/summary.txt
	echo "> Time Spent: $(( ( $TIME_END - $TIME_START ) / 1000 )) seconds" | tee -a outputs/temp/summary.txt
	echo "#" >> outputs/temp/summary.txt
#else
	echo "> UNABLE TO START BASELINE SYMMETRY TEST" >> outputs/temp/summary.txt
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