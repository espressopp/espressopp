#!/bin/bash - 
#===============================================================================
#
#          FILE: run_diff.sh
# 
#         USAGE: ./run_diff.sh 
# 
#   DESCRIPTION: 
# 
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: Federico G. Padua (), padua@uni-mainz.de
#  ORGANIZATION: Johannes Gutenberg Universitat Mainz, ZDV
#       CREATED: 02/11/2016 19:04
#      REVISION:  ---
#===============================================================================

FILE1=$1
FILE2=$2
FILE3=$3

OUT1=$(diff $1 $3)
echo "Should expect no output from diff, since we match the current time step"
if [ -z "$OUT1" ]; then
    echo "OUT1 is empty"
fi
echo "$OUT1"


OUT2=$(diff $2 $3)
echo "Should expect the second line (comment line) to differ, since we match the current time step"
if [ -n "$OUT2" ]; then
    echo "OUT2 not empty"
fi
echo "$OUT2"

