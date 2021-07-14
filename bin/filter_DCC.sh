#!/usr/bin/bash

input=$1
base=$(basename $input .txt)

## filter low reads
awk '{if($5 > 1) print $0}' $input > ${base}.filtered

## fix start position (+1) compared to circexplorer2, find_circ, circRNA_finder
awk -v OFS="\t" '{$2-=1;print}' ${base}.filtered > ${base}.bed
