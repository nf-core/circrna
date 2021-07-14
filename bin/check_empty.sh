#!/usr/bin/bash

## When consolidating tool outputs for a sample, check if any are empty.
## This will prevent R from throwing errors

for i in *bed; do

	if [[ -s $i ]];
	then
		:
	else
		echo "${i}" >> remove_empty_sample.txt
	fi

done

## sample.csv contains all tool bed files by default 
mv samples.csv checkme.csv

## if there are empty tool outputs then remove them from samples.csv and rewrite.
if [[ -s remove_empty_sample.txt ]];
then
	grep -vf remove_empty_sample.txt checkme.csv > samples.csv
else
	mv checkme.csv samples.csv
fi
