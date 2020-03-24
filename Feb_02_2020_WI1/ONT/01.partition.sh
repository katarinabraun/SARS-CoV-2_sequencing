#! /bin/bash

# demultiplex partitioned ONT FASTQ sequences

## step 1 - create gzip-compressed FASTQ file from all fastq_pass reads ##
# this will enable partitioning into an arbitrary number of FASTQ files
# the standard 4000 sequences per file contains too few reads to run efficiently on CHTC

FASTQ_PASS_FOLDER=$1

# compress files in FASTQ_PASS_FOLDER to single FASTQ file
mkdir -p merged_fastq
cat $FASTQ_PASS_FOLDER/*/*.fastq | pigz > merged_fastq/merged.fastq.gz

## step 2 - partition merged FASTQ into specified number of smaller FASTQ files ##
NUMBER_PARTITIONS=$2

mkdir -p partitioned_fastq

partition.sh \
qin=33 \
ow=t \
in=merged_fastq/merged.fastq.gz \
out=partitioned_fastq/%.fastq \
ways=$NUMBER_PARTITIONS

