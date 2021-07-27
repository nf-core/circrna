#!/usr/bin/bash 

input=$1
base=$(basename $input .gtf)

## strip headers from file and remove BSJ with less than 2 reads
grep -v "#" $input | grep -v "bsj 1.000" > ${base}.filtered
awk '{print $14}' ${base}.filtered | cut -d'.' -f1 > counts
cut -f 1,4,5,7 ${base}.filtered > ${base}.txt
paste ${base}.txt counts > ${base}.bed.tmp

## fix start position (+1) fix to match other tools
awk -v OFS="\t" '{$2-=1;print}' ${base}.bed.tmp > ${base}.bed
