#! /bin/bash

# merge demultiplexed samples to one FASTQ files per barcode
DEMULTIPLEXED_PATH='demultiplexed'

# make file for barcodes 1-24, even though many of these will be empty

mkdir -p merged_demultiplexed

for i in $(seq -f "%02g" 24 $END); 
    do 
        cat $DEMULTIPLEXED_PATH/*/barcode$i.fastq | pigz > merged_demultiplexed/barcode$i.fastq.gz
        echo $i;
    done

# remove empty files
for f in merged_demultiplexed/*
do
    if [[ $(gunzip -c $f | head -c1 | wc -c) == "0" ]] 
    then
        rm $f
    fi
done