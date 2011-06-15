#!/bin/bash
# Loops until all the S### directories contain the file FINISHED
# Then it exits

ndir=`ls -1 -d S[0-1][0-9][0-9] | wc -l`
dirs=`ls -d S[0-1][0-9][0-9]`

count=0

while [ $count -ne $ndir ]
do
   count=0
   for dir in $dirs
   do
      echo $dir
      if [ -e $dir/FINISHED ]
      then
         let count=count+1
      fi
   done
   sleep 5
done
