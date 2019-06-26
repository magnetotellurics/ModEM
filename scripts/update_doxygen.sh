#!/bin/bash
#
# NO ARGUMENTS: 
#
#
mkdir -p docs/
#
# STRING NOW
now=$(date "+%Y/%m/%d - %H:%M:%S")
#
#
echo "### START UPDATE DOXYGEN AT $now ###" >> outputs/update_oxygen.txt
#
#
cd src/
#
#
doxygen ../docs/doxygen_modem_on_config &>> ../outputs/update_oxygen.txt
#
#
echo "### FINISH DOXYGEN ###" >> ../outputs/update_oxygen.txt
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
# END OF SCRIPT

