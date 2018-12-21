#!/bin/bash

### DESCRIPTION: Rename indices of file names from library_ID to experiment_ID

dir='~/output/bdg/'

refile='~/sources/convert_libID_to_exp_ID.txt'

cd $dir
while read line; do
  file=$(echo $line | awk -F " " '{print $1}')
  index=$(echo $line | awk -F " " '{print $2}')
  finder=$(find -name ${file}*)
  if [ -f $finder ]; then
    idx=$(echo $finder | awk -F "_" '{print $1}')
    mv $finder ${index}${finder#*${idx}}
  fi
  
done < $refile
