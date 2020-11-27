#!/usr/bin/bash

file=$1

awk '{print $11}' $file | tr ',' '\n' | awk 'NR>1{print buf}{buf=$0}' > exon_sizes
seq $(awk '{print $10}' $file) > chrom
paste chrom exon_sizes > tmp.txt

rm -f exon_sizes
rm -f chrom

awk '$1="exon"$1' tmp.txt > tmp2.txt
awk -v OFS="\t" '{print $1, 0, $2}' tmp2.txt > circlize_exons.txt

rm -f tmp.txt
rm -f tmp2.txt
