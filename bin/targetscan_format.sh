#!/usr/bin/env

## Script that converts miRbase (mature.fa) file to
## TargetScan compatability. The motivation for doing
## this is the mature.fa file contains many more
## species than TargetScans miR_Family_Info.txt file.

## Stragtegy is simply to output a 3 column tab delim
## text file containing:
## 1. miR ID
## 2. miR (7bp) seed sequence from mature seq
## 3. Species ID (set to 0000, not important for output).

## Subset mature.fa according to the species provided by user to '--genome'.

## Stage input mature.fa file, species
MATURE="$1"
GENOME_ID="$2"

## Uncompress if necessary
if [ ${MATURE: -3} == ".gz" ]; then
    gunzip -f $MATURE
    MATURE=${MATURE%%.gz}
fi

## Convert to TargetScan
## Isolate the sequences
grep -v ">" $MATURE > mature_sequence
## Extract seed sequence (7bp after 1st)
awk '{print substr($1, 2, 7)}' mature_sequence > seed_sequence
## Isolate ID (awk last field (NF))
grep ">" $MATURE | awk -F ' ' '{print $NF}' > miR_ID
## Combine
paste miR_ID seed_sequence > targetscan_tmp.txt
## Correct delimiter, add dummy species
awk -v OFS="\t" '{print $1, $2, "0000"}' targetscan_tmp.txt > ${GENOME_ID}_targetscan.txt

## Tidy the work dir (uncomment these for debugging scratch dirs)
rm -rf mature_sequence
rm -rf miR_ID
rm -rf targetscan_tmp.txt
rm -rf seed_sequence
