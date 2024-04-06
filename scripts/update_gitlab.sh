#!/bin/bash
#
# ARGUMENTS: 1 - BRANCH NAME
#
#
branch_name=$1
# STRING NOW
now=$(date "+%Y/%m/%d - %H:%M:%S")
#
#
echo "#### START UPDATE $branch_name FROM AT $now ####"
#
#
git config --global user.name "protew"
git config --global user.email "paulowerdt@gmail.com"
git add .
git commit -m "GitLab '$branch_name' Runner Push [skip ci]"
git push https://protew:BGgwESV8qpGsBBdUZ8yk@gitlab.com/on.multiphysics/modem-on.git HEAD:$branch_name
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
echo "#### FINISH UPDATE GITLAB FROM AT $now ####"
#
#
exit 0
#
# END OF SCRIPT