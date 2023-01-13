#
# JMult_test 4 TXS
#
mpirun -np 9 ./Mod3DMT_STD -F ../../modem-oo/inputs/JMult_test/pr.ws ../../modem-oo/inputs/JMult_test/d.txt ON_FWD_JMULT_TEST_PRED_DATA.txt wEsoln1 ../../modem-oo/inputs/Others/modem_on_control_file.txt
#
mpirun -np 9 ./Mod3DMT_STD -T ../../modem-oo/inputs/JMult_test/pr.ws ../../modem-oo/inputs/JMult_test/d.txt ON_JT_JMULT_TEST_DSIGMA.txt ../../modem-oo/inputs/Others/modem_on_control_file.txt
#
mpirun -np 9 ./Mod3DMT_STD -M ../../modem-oo/inputs/JMult_test/pr.ws ../../modem-oo/inputs/JMult_test/dm.ws ../../modem-oo/inputs/JMult_test/d.txt ON_J_JMULT_TEST_JMHAT.txt ../../modem-oo/inputs/Others/modem_on_control_file.txt
#
mpirun -np 9 ./Mod3DMT_STD -I DCG ../../modem-oo/inputs/JMult_test/pr.ws ../../modem-oo/inputs/JMult_test/d.txt ../../modem-oo/inputs/Others/modem_on_control_file.txt ../../modem-oo/inputs/Others/modem_on_control_file.txt
#
# 1st_Example 16 TXS
#
mpirun -np 33 ./Mod3DMT_STD -F  ../../modem-oo/inputs/1st_Example/rFile_Model ../../modem-oo/inputs/1st_Example/rFile_Data_MT ON_FWD_1ST_EXAMPLE_PRED_DATA.txt wEsoln2 ../../modem-oo/inputs/Others/modem_on_control_file.txt
#
mpirun -np 33 ./Mod3DMT_STD -T  ../../modem-oo/inputs/1st_Example/rFile_Model ../../modem-oo/inputs/1st_Example/rFile_Data_MT ON_JT_1ST_EXAMPLE_DSIGMA.txt ../../modem-oo/inputs/Others/modem_on_control_file.txt
#
mpirun -np 33 ./Mod3DMT_STD -M ../../modem-oo/inputs/1st_Example/rFile_Model ../../modem-oo/inputs/1st_Example/rFile_Model ../../modem-oo/inputs/1st_Example/rFile_Data_MT ON_J_1ST_EXAMPLE_JMHAT.txt ../../modem-oo/inputs/Others/modem_on_control_file.txt
#
mpirun -np 33 ./Mod3DMT_STD -I DCG ../../modem-oo/inputs/1st_Example/rFile_Model ../../modem-oo/inputs/1st_Example/rFile_Data_MT ../../modem-oo/inputs/Others/modem_on_control_file.txt ../../modem-oo/inputs/Others/modem_on_control_file.txt