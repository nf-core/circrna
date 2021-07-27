#!/usr/bin/bash

## Stage files 

fasta=$1
input_file=$2

## Remove circRNA with BSJ reads < 1

awk '{if($13 > 1) print $0}' $input_file > filtered.circrna.txt

## Generate circRNA fasta file 

base=$(echo $input_file | cut -f 1 -d.)
counter=0

        while read line
        do
                chr=$(echo $line | awk '{ print $1 }')
                start=$(echo $line | awk '{ print $2 }')
                stop=$(echo $line | awk '{ print $3 }')
                strand=$(echo $line | awk '{ print $6 }')
                sam_out=$(samtools faidx $fasta "$chr:$start-$stop")
                header=">${chr}:${start}-${stop}_${strand}"
                echo $header >> ${base}.fa
                if [ $strand == "+" ]; then
                        sequence=$(echo "$sam_out" 2>&1 | tail -n +2 )
                elif [ $strand == "-" ]; then
                        sequence=$(echo "$sam_out" 2>&1 | tail -n +2 | rev | tr {AGTCagtc} {TCAGtcag})
                fi
                        echo "$sequence" >> ${base}.fa
                counter=$(($counter+1))
                if [ $(($counter%100)) -eq 0 ]; then
                        echo $counter
                fi

        done < filtered.circrna.txt
