#!/usr/bin/env bash

file=$1

while IFS='' read -r line; do

    name=$(echo $line)
    chr=$(echo $line | cut -d: -f1)
    start=$(echo $line | cut -d- -f1 | cut -d: -f2)
    stop=$(echo $line | cut -d- -f2 | cut -d: -f1)
    sign=$(echo $line | cut -d: -f3)

    echo -e "$chr\t$start\t$stop\t$name\t0\t$sign" >> ${name}.bed

done < $file
