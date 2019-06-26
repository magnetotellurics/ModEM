#!/bin/bash
#
# ARGUMENTS: 1 - JOB NAME
job_name=$1
#
# STRING NOW
now=$(date "+%Y/%m/%d - %H:%M:%S")
#
#
echo "### START UPDATE GITLAB FROM '${job_name}' AT $now ###" >> outputs/update_gitlab.txt
#
#
git config --global user.name "protew" &>> outputs/update_gitlab.txt
git config --global user.email "paulowerdt@gmail.com" &>> outputs/update_gitlab.txt
git add . &>> outputs/update_gitlab.txt
git commit -m "GitLab Runner Push [skip ci]" &>> outputs/update_gitlab.txt
git push https://protew:BGgwESV8qpGsBBdUZ8yk@gitlab.com/on.multiphysics/modem-on.git HEAD:master &>> outputs/update_gitlab.txt
#
#
# END OF SCRIPT
