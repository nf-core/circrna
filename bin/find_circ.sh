#!/bin/bash

# nextflow usage
# find_circ/find_circ.sh $genome $index <prefix.*.bwt2> $fastq_baseName $fastq[0] $fastq[1]

genome=$1
index=$2
base=$3
R1=$4
R2=$5



        bowtie2 -p 8 --very-sensitive --mm -D 20 --score-min=C,15,0 \
        -x $index -q -1 $R1 -2 $R2 \
        | samtools view -hbuS - | samtools sort --threads 8 -m 2G - > ${base}.bam

        samtools view -hf 4 ${base}.bam | samtools view -Sb - > ${base}_unmapped.bam

        unmapped2anchors.py ${base}_unmapped.bam | gzip > ${base}_anchors.qfa.gz

        bowtie2 -p 8 --reorder --mm -D 20 --score-min=C,-15,0 -q -x $index \
        -U ${base}_anchors.qfa.gz | find_circ.py -G $genome -p ${base} -s ${base}.sites.log > ${base}.sites.bed 2> ${base}.sites.reads

        echo "#chrom:start:end:name:n_reads:strand:n_uniq:best_qual_A:best_qual_B:spliced_at_begin:spliced_at_end:tissues:tiss_counts:edits:anchor_overlap:breakpoints" > $base/tmp

        cat $base/tmp | tr ':' '\t' > ${base}.circ_candidates.bed

        rm $base/tmp

        grep circ ${base}.sites.bed | grep -v chrM | sum.py -2,3 | scorethresh.py -16 1 | scorethresh.py -15 2 | scorethresh.py -14 2 | scorethresh.py 7 2 | scorethresh.py 8,9 35 | scorethresh.py -17 100000 >> ${base}.circ_candidates.bed

        rm ${base}.bam
        rm ${base}_anchors_unmapped.bam
        rm ${base}_anchors.qfa.gz
