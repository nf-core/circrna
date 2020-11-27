#!/usr/bin/bash

mkdir -p parent_genes

while IFS='' read -r line; do
        name=$(echo $line | awk '{print $4}')
        touch ${name}.bed
        echo "$line" >> ${name}.bed_tmp
        sed 's/[\t]*$//' ${name}.bed_tmp > ${name}.bed && rm ${name}.bed_tmp
        gene=$(bedtools intersect -a ${name}.bed -b filt.gtf -s -f 1.00 -wb | awk -F'gene_name ' '{print $2}' | awk -F';' '{print $1}' | sed 's/"//g' | uniq | head -n 1)
        rm ${name}.bed

        echo -e "$name\t$gene" >> parent_genes/${name}.parent_genes.txt

done < de_circ.bed
