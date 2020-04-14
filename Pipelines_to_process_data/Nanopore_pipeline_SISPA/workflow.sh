#! /bin/bash

# run through complete workflow using a combination of shell scripts and snakemake files

## Parameters ##
FASTQ_PASS_FOLDER='fastq_pass'
NUMBER_FASTQ_PARTITIONS='36'
MINIMUM_READ_LENGTH='300'


## partition sequences ##
bash /01.partition.sh  $FASTQ_PASS_FOLDER $NUMBER_FASTQ_PARTITIONS


## demultiplex sequences ##
snakemake \
--snakefile /02.demultiplex.snakefile \
--config \
min_length=$MINIMUM_READ_LENGTH \
partitioned_fastq_folder=partitioned_fastq \
--cores $NUMBER_FASTQ_PARTITIONS


## merge demultiplexed sequences into one FASTQ per barcode ##
bash /03.merge-demultiplex.sh \


## Trim out barcode sequences and get rid of LQ reads
snakemake --snakefile /04.subsample_QC.snakefile --config merged_demultiplexed=merged_demultiplexed --cores 12


## remove host and reagent reads ##
# this step is run separately on each barcode because the host databases may be different when multiple samples from multiple species are run in a single ONT run
preprocess () {
    snakemake \
    --snakefile /05.remove-host-reagent.snakefile \
    --config \
    ont_fastq_gz=$1 \
    reagent_db=/22592-reagent-db.fasta.gz \
    host_rna_db=$2 \
    host_dna_db=$3 
}

preprocess subsample/barcode##.fastq.gz PATH_TO_HOST_RNA_FILE  PATH_TO_HOST_RDNA_FILE;
preprocess subsample/barcode##.fastq.gz PATH_TO_HOST_RNA_FILE  PATH_TO_HOST_RDNA_FILE;
preprocess subsample/barcode##.fastq.gz PATH_TO_HOST_RNA_FILE  PATH_TO_HOST_RDNA_FILE;
preprocess subsample/barcode##.fastq.gz PATH_TO_HOST_RNA_FILE  PATH_TO_HOST_RDNA_FILE;


## Map cleaned reads to reference genome 
## This step is run separately on each barcode because the reference genome may be different between samples 
preprocess () {
    snakemake \
    --snakefile /06.map-reference-genome.snakefile \
    --config \
    ont_fastq_gz=$1 \
    mapping_genome=$2
}

preprocess cleaned/barcode##.clean.fastq.gz PATH_TO_REFERENCE_FULL_GENOME_FASTA
preprocess cleaned/barcode##.clean.fastq.gz PATH_TO_REFERENCE_FULL_GENOME_FASTA
preprocess cleaned/barcode##.clean.fastq.gz PATH_TO_REFERENCE_FULL_GENOME_FASTA
preprocess cleaned/barcode##.clean.fastq.gz PATH_TO_REFERENCE_FULL_GENOME_FASTA


## Converts mapped bam files to fastq files of only mapped reads
snakemake --snakefile /06.bam_to_fastq.snakefile --config mapped=mapped mapping_genome=PATH_TO_REFERENCE_FULL_GENOME_FASTA --cores 12


## Map cleaned reads to reference genome broken up into per genes 
preprocess () {
    snakemake \
    --snakefile /07.map-by-gene.snakefile \
    --config \
    ont_fastq_gz=$1 \
    mapping_genome=$2
}

preprocess cleaned/barcode##.clean.fastq.gz PATH_TO_REFERENCE_GENOME_BY_GENE_FASTA
preprocess cleaned/barcode##.clean.fastq.gz PATH_TO_REFERENCE_GENOME_BY_GENE_FASTA
preprocess cleaned/barcode##.clean.fastq.gz PATH_TO_REFERENCE_GENOME_BY_GENE_FASTA
preprocess cleaned/barcode##.clean.fastq.gz PATH_TO_REFERENCE_GENOME_BY_GENE_FASTA
preprocess cleaned/barcode##.clean.fastq.gz PATH_TO_REFERENCE_GENOME_BY_GENE_FASTA


## Calls variants by gene using a fastq that is separated by the gene name 
preprocess () {
    snakemake \
    --snakefile /08.call-variants-by-gene.snakefile \
    --config \
    bygenebam=$1 \
    mapping_genome=$2
}

preprocess mapped_bygene/barcode##.primary.bam PATH_TO_REFERENCE_GENOME_BY_GENE_FASTA
preprocess mapped_bygene/barcode##.primary.bam PATH_TO_REFERENCE_GENOME_BY_GENE_FASTA
preprocess mapped_bygene/barcode##.primary.bam PATH_TO_REFERENCE_GENOME_BY_GENE_FASTA
preprocess mapped_bygene/barcode##.primary.bam PATH_TO_REFERENCE_GENOME_BY_GENE_FASTA


## make minnashes from cleaned reads vs. nt database ##
snakemake --snakefile /07.minhash.dataset.snakefile --config cleaned=cleaned --cores 12
snakemake --snakefile /07.minhash.sequences.snakefile --config cleaned=cleaned --cores 12


## cleanup
# rm -rf demultiplexed merged_demultiplexed merged_fastq partitioned_fastq tmp