#!/bin/bash
#
# NO ARGUMENTS: 
#
echo "#### START UPDATE DOXYGEN DOCUMENTATION FOR ModEM-OO AT $(date "+%Y/%m/%d - %H:%M:%S")" | tee -a outputs/temp/summary.txt
echo "   #" | tee -a outputs/temp/summary.txt
#
cd src/
#
#
doxygen ../docs/doxygen_modem_oo_config
#
# CATCH RESULT
result=$?
#
# TEST RESULT
if [ "$result" -ne "0" ]; then
	#
	#
	echo "update_oxygen FAIL: $result"
	#
	#
	cd ..
	#
	#
	exit $result
fi
#
#
echo "#### FINISH UPDATE DOXYGEN ####"
#
# REMOVE OLD html/ AT docs/
rm -rf ../docs/html
#
# SEND NEW html/ TO docs/
mv html/ ../docs/
#
# REMOVE NEW latex/ CREATED
rm -rf latex/
#
#
cd ..
#
# 
exit 0
#
# END OF SCRIPT

