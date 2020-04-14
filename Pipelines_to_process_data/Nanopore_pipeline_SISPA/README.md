# GoingViral
GoingViral is a pipeline for the analysis of SISPA generated viral metagonomics reads. It's a configurable pipeline that performs multiple functions to generate clean viral metagenomics reads that can be used for downstream analysis. It requires minimal dependencies as the pipeline relies on a docker container to host the software.

## Table of Contents
* [Requirements](#requirements)
* [Installing](#installing)
* [Contents](#contents)
* [Usage](#usage)
* [Running Workflow](#running)

### Requirements
* Linux or MacOS
* Docker


## Installing
* Install the [Docker CE engine](https://docs.docker.com/install/)

1. Change to directory of fastq_pass files from ONT sequencer:
`
cd PATH_TO_FILES
`
2. Pull the most up to date docker container:
    ```docker pull gkmoreno/goingviral:v1```
3. Launch the docker container
    ```docker run -it -v $(pwd):/scratch -w /scratch gkmoreno/goingviral:v1 /bin/bash```


### Contents
The entire GoingViral workflow is uploaded as a series of snakemake and bash scripts that can be run one by one or can by run sequentially using the `workflow.sh` script as a driver script. All of these bash and snakemake scripts  except for the `workflow.sh` will be provided in the docker container. 
A brief description of the GoingViral Docker Container contents:
- `01.partition.sh` - will combine all ONT fastq_pass files into one merged folder and then will partition it out into 36 sub-folders to demultiplex simultaenously. 
- `02.demultiplex.snakefile` - Runs `qcat` on the 36 sub-folders. Discards reads <300bp in length. Trims out ONT adaptors and barcodes. 
- `03.merge-demultiplex.sh` - Merges the demultiplexed reads into a single fastq.gz for each barcode using `pigz`.
- `04.subsample_QC.snakefile` - Discards reads ≤Q7 and trims out SISPA primer sequence using `reformat.sh`.
- `05.remove-host-reagent.snakefile` - Uses `minimap` to bioinformatically deplete of host and reagent contaminants.
- `06.map-reference-genome.snakefile` - Maps cleaned reads to a reference file using `minimap` and will call variants ≥10% frequency using `callvariants.sh`. 
- `06.bam_to_fastq.snakefile` - Converts the mapped bam file to fastq using `reformat.sh`.
- `07.map-by-gene.snakefile` - Maps cleaned reads to a reference file composed of only the coded gene regions using `minimap`. 
- `08.call-variants-by-gene.snakefile`- Call variants in the coded gene regions using `callvariants.sh`. 
- `09.minhash.dataset.snakefile` - Performs `sendsketch.sh` on the entire dataset - outputs the top 1000 hits in your sample. 
- `09.minhash.sequences.snakefile` - Performs `sendsketch.sh` and classifies reads on a per sequence basis 

### Usage
GoingViral uses a bash script (workflow.sh) to provide parameters to the pipeline. To run the pipeline simply upodate the workflow.sh script to use the configurations that you want.
Areas to configure: 
* PATH_TO_HOST_RNA_FILE
* PATH_TO_HOST_DNA_FILE
* PATH_TO_REFERENCE_FULL_GENOME_FASTA
* PATH_TO_REFERENCE_GENOME_BY_GENE_FASTA

```
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
```

### Running
Once your workflow.sh has been configured and docker container has been launched - you can start the workflow by simply running:
`bash worklfow.sh`
in the terminal window. 
