#!/usr/bin/env bash

## .fa = backsplice .fasta = canonical to publish

file_prefix=$(basename $1 .fa)

# Split multi-fasta file into single records
mkdir canonical_seqs && mv $1 canonical_seqs/
ls -d canonical_seqs/*.fa | while read -r line; do
  cat $line | awk '{if (substr($0, 1, 1)==">") {filename=(substr($0,2) ".fa")} print $0 > filename}'
done

# add backsplice site to individual fasta files
for file in *.fa; do
  fn=$(basename $file .fa)
  cat $file | grep ">" > header.txt
  cat $file | grep -v ">" | cut -c1-20 > 20bp.txt
  cat $file | grep -v ">" > seqs.txt
  paste seqs.txt 20bp.txt | sed -e 's/\t//' > seqs_20bp.txt
  paste header.txt seqs_20bp.txt | sed -e 's/\t/\n/' > ${fn}.fasta
  rm header.txt 20bp.txt seqs.txt seqs_20bp.txt
done

# rm single fasta entries from while read line
rm *.fa
# merge backsplice fastas into one for miRNA pred
cat *.fasta > ${file_prefix}.fa
# remove intermediate for loop fasta file
rm *.fasta
# move canonical seqs back to publish to outdir
ls -d canonical_seqs/*.fa | while read -r line; do
  mv $line ${file_prefix}.fasta
done
