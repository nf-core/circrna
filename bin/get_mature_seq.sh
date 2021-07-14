#!/usr/bin/bash


echo "===================================================================================="
echo "[nf-core/circrna]: circRNA annotation script                                        "
echo "[nf-core/circrna]: Author: Barry Digby                                              "
echo "[nf-core/circrna]: Institution: National University of Ireland, Galway              "
echo "===================================================================================="

mkdir -p bed12

while IFS='' read -r line; do

    name=$(echo $line | awk '{print $4}')
    touch ${name}.bed
    echo "$line" >> ${name}.bed_tmp
    sed 's/[\t]*$//' ${name}.bed_tmp > ${name}.bed && rm ${name}.bed_tmp

    bedtools intersect -a filt.gtf -b ${name}.bed -s -f 1.00 > ${name}.gtf

    start=$(echo $line | awk '{print $2}')
    stop=$(echo $line | awk '{print $3}')

    echo "[nf-core/circrna]: Starting analysis for: $name"

    # is the gtf file NOT (-s) empty? i.e did it overlap biotypes?
    if [[ -s ${name}.gtf ]];
    then

        echo "[nf-core/circrna]: $name overlaps features in GTF file"
        echo "[nf-core/circrna]: Inspecting Genes..."

        gene_id=$(awk -F'gene_id ' '{print $2}' ${name}.gtf | \
        awk -F';' '{print $1}' | sed 's/"//g' | uniq | paste -s -d, -)  ## remove uniq here?

        echo "[nf-core/circrna]: Overlapping Gene IDs: $gene_id"
        echo "[nf-core/circrna]: Converting to BED12"

        gtfToGenePred ${name}.gtf ${name}.genepred
        genePredToBed ${name}.genepred ${name}_predtobed.bed

        # Attempting perfect exon boundary overlaps
        echo "[nf-core/circrna]: Attempting to fit circRNA to exon boundaries"
        awk -v OFS="\t" -v start="$start" -v stop="$stop" \
        '{if($2==start && $3==stop) print $0}' ${name}_predtobed.bed | \
        sort -rnk10 | head -n 1 > ${name}.bed12.bed

        # Resulting file not empty? i.e perfectly overlapped with exon boundaries?
        if [[ -s ${name}.bed12.bed ]];
        then
            echo "[nf-core/circrna]: ${name} overlaps exons, is a circRNA"
            type="circRNA"
        else

            echo "[nf-core/circrna]: circRNA imperfectly overlaps exons"
            echo "[nf-core/circrna]: Investigating if EIciRNA or acceptable to take underlying transcript"
            echo "[nf-core/circrna]: Retrying with longest underlying transcript"

            awk -v OFS="\t" '{$13 = $3 - $2; print}' ${name}_predtobed.bed | \
            sort -rnk13 | cut -f13 --complement | head -n 1 > ${name}.bed12.bed_tmp

            tx_len=$(awk -v OFS="\t" '{$13 = $3 - $2; print}' ${name}_predtobed.bed | \
            sort -rnk13 | awk '{print $13}' | head -n 1)

            circ_len=$(awk -v OFS="\t" '{$7 = $3 - $2; print}' ${name}.bed | awk '{print $7}')

            echo "[nf-core/circrna]: Best transcript length: $tx_len"
            echo "[nf-core/circrna]: $name length: $circ_len"

            difference=$(($circ_len - $tx_len))

            if [[ $difference -gt 200 ]];
            then

                echo "[nf-core/circrna]: Transcript exon boundaries more than 200nt off $name"
                echo "[nf-core/circrna]: Treating as EIciRNA"

                type="EIciRNA"
                block_count=1
                block_size=$(($stop-$start))
                rgb="0,0,0"
                block_start=0
                awk -v OFS="\t" -v thick=$start -v rgb=$rgb -v count=$block_count -v start=$block_start -v size=$block_size \
                '{print $0, thick, thick, rgb, count, size, start}' ${name}.bed > ${name}.bed12.bed
                rm ${name}.bed12.bed_tmp
            else

                echo "[nf-core/circrna]: Transcript exon boundaries within 200nt of ${name}"
                echo "[nf-core/circrna]: Treating ${name} as circRNA."
                type="circRNA"
                mv ${name}.bed12.bed_tmp ${name}.bed12.bed
            fi
        fi
        else

            echo "[nf-core/circrna]: $name returned no GTF overlaps."
            echo "[nf-core/circrna]: Treating as an intronic circRNA"

            type="ciRNA"
            block_count=1
            block_size=$(($stop-$start))
            rgb="0,0,0"
            block_start=0
            awk -v OFS="\t" -v thick=$start -v rgb=$rgb -v count=$block_count -v start=$block_start -v size=$block_size \
            '{print $0, thick, thick, rgb, count, size, start}' ${name}.bed > ${name}.bed12.bed
    fi

    awk -v type="$type" 'BEGIN{FS=OFS="\t"}{$9=type}1' ${name}.bed12.bed > ${name}.bed12.bed_tmp
    awk -v OFS="\t" -v name=$name '{$4 = name; print}' ${name}.bed12.bed_tmp > ${name}.bed12.bed_tmp1

    rm ${name}.bed12.bed
    rm ${name}.bed12.bed_tmp
    mv ${name}.bed12.bed_tmp1 ${name}.bed12.bed
    echo "[nf-core/circrna]: cleaning up intermediate files"
    rm -f ${name}.gtf
    rm -f ${name}.genepred
    rm -f ${name}_predtobed.bed
    rm -f ${name}.bed

    cp ${name}.bed12.bed bed12/
    rm -rf ${name}.bed12.bed
    echo "===================================================================================="
    echo "===================================================================================="

done < circs.bed

cat bed12/*.bed12.bed > master_bed12.bed
