#!/usr/bin/bash

mkdir -p bed12

type

while IFS='' read -r line; do

    name=$(echo $line | awk '{print $4}')
    touch ${name}.bed
    echo "$line" >> ${name}.bed_tmp
    sed 's/[\t]*$//' ${name}.bed_tmp > ${name}.bed && rm ${name}.bed_tmp

    bedtools intersect -a filt.gtf -b ${name}.bed -s -f 1.00 > ${name}.gtf

    start=$(echo $line | awk '{print $2}')
    stop=$(echo $line | awk '{print $3}')

    echo "#"
    echo "## Starting analysis for: $name"
    echo "#"

    # is the gtf file NOT (-s) empty? i.e did it overlap biotypes?
    if [[ -s ${name}.gtf ]];
    then

        echo "$name overlaps features in GTF file"
        echo "Inspecting biotypes..."

        gene_types=$(awk -F'gene_type ' '{print $2}' ${name}.gtf | \
        awk -F';' '{print $1}' | sed 's/"//g' | uniq)

        echo "$gene_types"

        gtfToGenePred ${name}.gtf ${name}.genepred
        genePredToBed ${name}.genepred ${name}_predtobed.bed

        # Attempting perfect exon boundary overlaps
        awk -v OFS="\t" -v start="$start" -v stop="$stop" \
        '{if($2==start && $3==stop) print $0}' ${name}_predtobed.bed | \
        sort -rnk10 | head -n 1 > ${name}.bed12.bed

        # Resulting file not empty? i.e perfectly overlapped with exon boundaries?
        if [[ -s ${name}.bed12.bed ]];
        then
            type="circRNA"
        else

            echo "The circRNA imperfectly overlaps an exon"
            echo "Investigating if EIciRNA or acceptable to take longest transcript"
            echo "Retrying with longest transcript"

            awk -v OFS="\t" '{$13 = $3 - $2; print}' ${name}_predtobed.bed | \
            sort -rnk13 | cut -f13 --complement | head -n 1 > ${name}.bed12.bed_tmp

            echo "Checking best transcript with $name"

            tx_len=$(awk -v OFS="\t" '{$13 = $3 - $2; print}' ${name}_predtobed.bed | \
            sort -rnk13 | awk '{print $13}' | head -n 1)

            circ_len=$(awk -v OFS="\t" '{$7 = $3 - $2; print}' ${name}.bed | awk '{print $7}')

            echo "Best transcript length: $tx_len"
            echo "$name length: $circ_len"

            difference=$(($circ_len - $tx_len))

            if [[ $difference -gt 200 ]];
            then

                echo "Transcript is more than 200nt off $name"
                echo "Treating as EIciRNA"

                type="EIciRNA"
                block_count=1
                block_size=$(($stop-$start))
                rgb="0,0,0"
                block_start=0
                awk -v OFS="\t" -v thick=$start -v rgb=$rgb -v count=$block_count -v start=$block_start -v size=$block_size \
                '{print $0, thick, thick, rgb, count, size, start}' ${name}.bed > ${name}.bed12.bed
                rm ${name}.bed12.bed_tmp
            else

                echo "Transcript is within 200nt of ${name}"
                echo "Taking best transcript as coordinates"
                type="circRNA"
                mv ${name}.bed12.bed_tmp ${name}.bed12.bed
            fi
        fi
        else

            echo "$name returned empty GTF file in bedtools query."
            echo "Most likely an intronic circRNA"

            type="ciRNA"
            block_count=1
            block_size=$(($stop-$start))
            rgb="0,0,0"
            block_start=0
            awk -v OFS="\t" -v thick=$start -v rgb=$rgb -v count=$block_count -v start=$block_start -v size=$block_size \
            '{print $0, thick, thick, rgb, count, size, start}' ${name}.bed > ${name}.bed12.bed
    fi

    echo "replacing tx with circRNA & RGB with circRNA type"
    awk -v type="$type" 'BEGIN{FS=OFS="\t"}{$9=type}1' ${name}.bed12.bed > ${name}.bed12.bed_tmp
    awk -v OFS="\t" -v name=$name '{$4 = name; print}' ${name}.bed12.bed_tmp > ${name}.bed12.bed_tmp1

    rm ${name}.bed12.bed
    rm ${name}.bed12.bed_tmp
    mv ${name}.bed12.bed_tmp1 ${name}.bed12.bed
    echo "cleaning up intermediate files"
    rm -f ${name}.gtf
    rm -f ${name}.genepred
    rm -f ${name}_predtobed.bed
    rm -f ${name}.bed

    cp ${name}.bed12.bed bed12/

done < circs.bed

cat bed12/*.bed12.bed > de_circ_exon_annotated.bed
