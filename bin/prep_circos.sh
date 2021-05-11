#!/usr/bin/bash

file=$1

## capture the exon sizes column, sep by comma if multiple
awk '{print $11}' $file | tr ',' '\n' > exon_sizes_tmp
## check number of exons by counting n lines
cat exon_sizes_tmp
len=$(wc -l exon_sizes_tmp | awk '{print $1}')
## print buf removes empty last line inserted in process.
## be careful here, do not execute this if only one line present (EIciRNA, ciRNA)
#[[ $len -gt 1 ]] && awk 'NR>1{print buf}{buf=$0}' exon_sizes_tmp > exon_sizes || mv exon_sizes_tmp exon_sizes
if [[ $len -gt 1 ]] ; then awk 'NR>1{print buf}{buf=$0}' exon_sizes_tmp > exon_sizes; else mv exon_sizes_tmp exon_sizes; fi
## seq converts '2' to '1', '2'.. but does not work on just '1' (EIciRNA, ciRNA)
n_exons=$(awk '{print $10}' $file)
if [[ $n_exons -gt 1 ]]; then seq $n_exons > chrom; else printf $n_exons > chrom; fi
#[[ $n_exons -gt 1]] && seq $n_exons > chrom || printf $n_exons > chrom
paste chrom exon_sizes > tmp.txt

rm -f exon_sizes
rm -f chrom

awk '$1="exon"$1' tmp.txt > tmp2.txt
awk -v OFS="\t" '{print $1, 0, $2}' tmp2.txt > circlize_exons.txt

rm -f tmp.txt
rm -f tmp2.txt
