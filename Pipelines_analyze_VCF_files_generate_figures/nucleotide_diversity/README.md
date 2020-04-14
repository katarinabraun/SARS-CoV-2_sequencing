# Intragene diversity 
  
**Author**: Katarina Braun 

## Overview
- The purpose of this script is to generate genewise nucleotide diversity (π) using the SNPGenie script and then to gnerate plots to visualize and compare genewise πS (synonymous diversity) and πN (nonsynonymous diversity). 

## Dependencies 
- pandas
- numpy
- matplotlib
- glob
- os
- functools
- pandas
- itertools
- glob
- os
- sklearn 
- pysam
- random 
- seaborn 

First, this notebook will go over syntax to run (SNPGenie)[https://github.com/chasewnelson/SNPGenie]. SNPGenie is a collection of Perl scripts for estimating πN/πS, dN/dS, and gene diversity from next-generation sequencing (NGS) single-nucleotide polymorphism (SNP) variant data. 

WITHIN-POOL ANALYSIS. Use snpgenie.pl, the original SNPGenie. Analyzes within-sample πN/πS from pooled NGS SNP data. SNP reports (VCF) must each correspond to a single population/pool, with variants called relative to one reference sequence (one sequence in one FASTA file).

Here's the syntax for running SNPGenie that I am using: 

```bash
perl snpgenie.pl --vcfformat=4 --snpreport=path/to/SNPREPORT.vcf --fastafile=path/to/ref.fasta --gtffile=/path/to/GTF.gtf

```

SNPGenie will not take in a reference and GTF divided into genes, but we used a reference that is divided into gene segments to generate original SNV calls. 

- **ORF1a**: +265    
- **ORF1b**: +13,468  
- **S**: +21,559  
- **ORF3a**: +25,389  
- **E**: +26,241  
- **M**: +26,519  
- **ORF6**: +27,189  
- **ORF7a**: +27,390  
- **ORF8**: +27,890  
- **N**: +28,780  
- **ORF10**: +29,554

Then, this notebook will take raw diveristy (π) data and generate plots to visualize and compare genewise πS (synonymous diveristy) and πN(nonsynonymous diversity). 
- will plot π, pulled from `population_summary.txt` files -- fifth column 
- plot each gene within each sample

## Input: 

- `SNPGenie_ouput/sample/product_results.txt`



## Output:  

- Nucleotide diversity figure.