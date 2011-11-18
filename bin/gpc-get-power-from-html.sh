#!/bin/sh
#echo $#
if test $#=0
then 
cat | sed -f  html_split.sedcmds $1 | grep navy | grep -v green | sed -f html_clean.sedcmds
else
sed -f  html_split.sedcmds $1 | grep navy | grep -v green | sed -f html_clean.sedcmds
fi




