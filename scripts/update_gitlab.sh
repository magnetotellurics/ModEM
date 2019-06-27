#!/bin/bash
#
# NO ARGUMENTS
#
# STRING NOW
now=$(date "+%Y/%m/%d - %H:%M:%S")
#
#
echo "#### START UPDATE GITLAB FROM AT $now ####"
#
#
git config --global user.name "protew"
git config --global user.email "paulowerdt@gmail.com"
git add .
git commit -m "GitLab Runner Push [skip ci]"
git push https://protew:BGgwESV8qpGsBBdUZ8yk@gitlab.com/on.multiphysics/modem-on.git HEAD:master
#
# CATCH RESULT
result=$?
#
# TEST RESULT
if [ "$result" -ne "0" ]; then
	#
	#
	echo "update_gitlab FAIL: $result"
	#
	#
	exit $result
fi
#
#
exit 0
#
# END OF SCRIPT
