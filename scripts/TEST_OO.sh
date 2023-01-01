#
# JMult_test 4 TXS
#
./TestSerial -f -m ../inputs/JMult_test/pr.ws -d ../inputs/JMult_test/d.txt -pd OO_FWD_JMULT_TEST_PRED_DATA.txt -c ../docs/control_file_template
#
./TestSerial -jt -m ../inputs/JMult_test/pr.ws -d ../inputs/JMult_test/d.txt -pd OO_JT_JMULT_TEST_PRED_DATA.txt -dm OO_JT_JMULT_TEST_DSIGMA.txt -c ../docs/control_file_template
#
./TestSerial -j -m ../inputs/JMult_test/pr.ws -pm ../inputs/JMult_test/dm.ws -d ../inputs/JMult_test/d.txt -pd OO_J_JMULT_TEST_PRED_DATA.txt -gd OO_J_JMULT_TEST_JMHAT.txt -c ../docs/control_file_template
#
./TestSerial -i -m ../inputs/JMult_test/pr.ws -d ../inputs/JMult_test/d.txt -pd OO_INV_JMULT_TEST_PRED_DATA.txt -dm OO_INV_JMULT_TEST_DSIGMA.txt -c ../docs/control_file_template
#
# 1st_Example 16 TXS
#
./TestSerial -f -m ../inputs/1st_Example/rFile_Model -d ../inputs/1st_Example/rFile_Data_MT -pd OO_FWD_1ST_EXAMPLE_PRED_DATA.txt -c ../docs/control_file_template
#
./TestSerial -jt -m ../inputs/1st_Example/rFile_Model -d ../inputs/1st_Example/rFile_Data_MT -pd OO_JT_1ST_EXAMPLE_PRED_DATA.txt -dm OO_JT_1ST_EXAMPLE_DSIGMA.txt -c ../docs/control_file_template
#
./TestSerial -j -m ../inputs/1st_Example/rFile_Model -pm ../inputs/1st_Example/rFile_Model -d ../inputs/1st_Example/rFile_Data_MT -pd OO_J_1ST_EXAMPLE_PRED_DATA.txt -gd OO_J_1ST_EXAMPLE_JMHAT.txt -c ../docs/control_file_template
#
./TestSerial -i -m ../inputs/1st_Example/rFile_Model -d ../inputs/1st_Example/rFile_Data_MT -pd OO_INV_1ST_EXAMPLE_PRED_DATA.txt -dm OO_INV_1ST_EXAMPLE_DSIGMA.txt -c ../docs/control_file_template
#
