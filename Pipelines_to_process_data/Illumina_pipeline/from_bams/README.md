# Sniffles

Sniffles was written by Dr. Kelsey Florek and Joe Lalli and was modified for use here by Katarina Braun. 
Instructions for use have been modified here to ensure data, analyses, and figures can be reproduced.  

**This is the first script that should be run to replicate analyses and figure-generation for the Illumina data.**

`Sniffles` is a pipeline originally written for the analysis of influenza genomes, but we've modified it here for the analysis of SARS-CoV-2 viruses. It's a configurable pipeline that performs multiple functions to generate variant and consensus level information. It requires minimal dependencies as the pipeline relies on a docker container to host the software.

This file contains the config parameters as well as the commands necessary to replicate the analysis we've done for this paper.

**Importantly, this directory contains a `sniffles.py` script which takes in the bam files available in the `data_raw` directory.** The `-B` flag should be used to indicate that the input is alignment `bams`, instead of paired-end FASTQs. 


## Table of Contents
* [Requirements](#requirements)
* [Installing](#installing)
* [Usage](#usage)

### Requirements
* MacOS
* Python 3.6 or later
* Docker (update Docker before running `Sniffles` by typing `conda update Docker` into the terminal)

- numpy>=1.15.4
- pandas>=0.23.4
- Bio>=0.1.0
- docker_py>=1.10.6
- PyYAML>=5.1
- slackclient>=2.0.1
- datetime

### Installing
* Make sure you are using the correct version of python
* Use git to get the latest sniffles build `git clone https://github.com/k-florek/sniffles.git`
* Install the [Docker CE engine](https://docs.docker.com/install/)
* Install the required python libraries by pointing the pip installer to the requirements document in the sniffles project `pip3 install -r requirements.txt`

### Usage
Sniffles uses a configuration file to provide parameters to the pipeline. There is a default configuration file included in the sniffles project. To run the pipeline simple point the `sniffles.py` program at the correct config.yml and directory containing the raw Illumina reads.

```
 #####
#     #  #    #  #  ######  ######  #       ######   ####
#        ##   #  #  #       #       #       #       #
 #####   # #  #  #  #####   #####   #       #####    ####
      #  #  # #  #  #       #       #       #            #
#     #  #   ##  #  #       #       #       #       #    #
 #####   #    #  #  #       #       ######  ######   ####



usage: sniffles.py [-h] [-c config] [-i input] [-o output] [-t threads]

Pipeline to examine SNPs from raw illumina reads

optional arguments:
 -h, --help  show this help message and exit
 -c config   config file
 -i input    raw reads directory - defaults to working directory
 -o output   output directory - defaults to working directory
 -t threads  number of cpus to use for pipeline
 -B   indicates that input will be in the form of bam files 
```

### Here's how things should be set up to run this pipeline:

Before running this pipeline, navigate to the directory which contains 
- `sniffles.py`
- `config.yml`
- `[ref].fasta`
- `[gtf].gtf`
- `supporting_code/*`

Naming schemes for the input files are important for getting this pipeline to run. Specifically, the reference and GTF files should contain the same rootname as the rootname of the bam file that is being fed into the file. For example, when running the `primary_NP_swab` sample, the reference fasta should be named `primary_NP_swab.fasta` and the GTF should be named `primary_NP_swab.gtf`. 
- The ref and and gtf have to be updated with each sample. The file itself does not need to be updated, just the names of these files.
- These filenames should also be updated in the config file: 
  - `outdir: 'primary_NP_swab'`
  - `logfile: 'primary_NP_swab_logs.log'`
  - `referenceSequences: ["primary_NP_swab.fasta"]`
  - `gtfFileNames: ["primary_NP_swab.gtf"]`

Additionally, the bam file being fed into the pipeline, `primary_NP_swab.bam` should be contained within a directory that is named `primary_NP_swab` and the directory that contains this dictory should be fed into the command. For example, if the location of the bam is `Volumes/bams/primary_NP_swab/primary_NP_swab.bam`, the -i flag should be followed by `Volumes/bams`. 

The command should be run from the directory that contains `sniffles.py`, `config.yml`, etc. 

Only one file should be run at a time. The script errors out if you try to run more than one at a time. 

Example usage: 

```bash
python sniffles.py -c config.yml -B -i Volumes/bams
```

### Run the above command for each of the following files: 

For each of the following files: 
1. p1_vero76
2. p1_veroE6
3. p1_veroSTAT1KO
4. p2a_vero76
5. p2b_vero76
6. primary_NP_swab

## Input 

Raw fastQ files:  
`SARSCoV2_passage_MS/data_raw/Illumina_bams/*`

Config file (parameters outlined below):  
`SARS-CoV-2/SARSCoV2_passage_MS/data_pipelines/Illumina_pipeline/from_bams/config.yml`

Reference fasta:  
`SARSCoV2_passage_MS/data_pipelines/from_bams/ref.fasta`

GTF file:  
`SARSCoV2_passage_MS/data_pipelines/from_bams/ref.gtf`

## Output 

All output files are written to:  
`SARSCoV2_passage_MS/data_pipelines/Illumina_pipeline/from_bams/[sample_name]/*`

## Config file parameters  

The config file "sample", "ref" and "GTF" need to be updated to match sample name before running each sample. 

```Python
#configuration file for holding runtime configurations and parameters

# ---Pipeline Execution Parameters---
exec:
  #name of the output directory
  outdir: 'sample'
  #name of the logfile
  logfile: 'sample_logs.log'
  #list of reference sequences used for mapping 
  #(brackets[] indicate more than one reference. 
  #If used, inputs must be in a subfolder named for the reference sequence)
  referenceSequences: ["sample.fasta"]
  #minimum average depth after mapping
  minimumAverageDepth: 100
  #minimum percent of reference bases covered
  percentRefCovered: 100

  # ---replicate variables---
  #samples in replicate
  replicates: False
  #Notation for sample replicates: SampleX_R1.fasta 
  #eg. if your fastq is H1N1.1_R1, please write "Sample.1_R1"
  # or if fastq is H1N1Rep1_R1, please write "SampleRep1_R1"
  # or if fastq is r1_ZW456_R1, please write "r1_Sample_R1"
  replicateNotation: "Sample_rep1_R1"
  #variant calling using 'Varscan' or 'RePlow' or 'Compare' (only three options). 
  #'Compare' will run both and produce an excel file that places the results side by side. 
  #'Compare' does not work with more than one set of replicates.
  #***Compare is not currently supported.***
  callSNPs: 'Varscan'
  
  # ---pipeline functions---
  #filter out samples which don't meet minimum coverage threshold
  coverageFilter: False
  #consensus generation using VarScan
  generateConsensus: True
  #consensus generation for all samples of given reference using VarScan
  generatePopConsensus: True
  #map the reads to the consensus sequence using bowtie2
  mapToConsensus: False
  #annotate SNPs with Louise Moncla's custom annotator
  annotateSNPs: True
  #use SNPgenie to derive population statistics #we'll do this separately
  SNPgenier: False

  # ---data cleaning---
  #use unpaired reads, or only paired reads
  unpaired: False

  #coverage depth normalization using seqtk
  #normalize coverage across the sam or bam file by randomly downsampling to a set number of reads
  normalizeCoverage: False
  totalReads: 200000

# ---Quality Trimming Parameters---
trimmomatic:
  #remove adapter sequences from the reads
  removeAdapters: False
  #adapter file name included in trimmomatic
  adaptersFileName: "Nextera_XT_adapter.fa"
  #paired end reads
  paired: True
  #minmial read length
  minlength: 100
  #sliding window size for quality trim
  windowSize: 5
  #quality threashold for trim in sliding window
  qscore: 30

# ---Single Nucelotide Polymorhpism Calling Parameters---
snpcalling:
  #minimum coverage for SNP
  minCoverage: 100
  #minimum quality for SNP
  snpQualityThreshold: 30
  #minimum frequency for SNP
  snpFrequency: 0.01
  #minimum frequency for a consensus SNP in generating the consensus sequence
  consensusFrequency: 0.5

#Replow-specific settings
replow_settings:
  #minimum mapping quality
  mapquality: 30
  #influenza mutation rate
  mutrate:  4E-05
  
#The fundamental question being asked is "What is the pre-test probability of any given nucleotide being different than the reference sequence?"
#This number is derived in cancer from looking at the number of mutation per Mb in the tumor compared to somatic mutations. You can't really do that in viruses.
#The closest analogy would be looking at the number of SNPs present in one ferret after infection with a standard strain. I don't know if people have done this,
#Because mutation rate/replication cycle is so much easier to do.
#This paper more recently estimated 1.4x10^-5 mutations/site/48hrs in influenza A, which seems like a reasonable estimate for a course of influenza.
#source: Nobusawa E, Katsuhiki S. "Comparison of the Mutation Rates of Human Influenza A and B Viruses". Journal of Virology Mar 2006, 80 (7) 3675-3678; DOI: 10.1128/JVI.80.7.3675-3678.2006

# ---Post-SNP population analysis---
postprocessing:
  #location of gtf file
  gtfFileNames: ["sample.gtf"]
  #minimum frequency of SNPs to be considered by SNPGenie
  minSNPfreq: 0.01
  #sliding window size for SNPGenie analysis (# of codons)
  slidingwindow: 50
```

