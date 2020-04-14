FROM continuumio/miniconda3

# update conda 
RUN conda upgrade -n base conda

#setup channels 
RUN export PATH="/opt/conda/bin:${PATH}" && \
conda config --add channels defaults && \
conda config --add channels bioconda && \ 
conda config --add channels conda-forge

# install tools from conda 
RUN conda install -y \
seqtk qcat minimap2 canu seqkit samtools bedtools pigz snakemake bbmap 

# ADD files from folder and push to dockerfile - no need to move files around 
ADD 01.partition.sh /01.partition.sh
ADD 02.demultiplex.snakefile /02.demultiplex.snakefile
ADD 03.merge-demultiplex.sh /03.merge-demultiplex.sh
ADD 04.subsample_QC.snakefile /04.subsample_QC.snakefile
ADD 05.remove-host-reagent.snakefile /05.remove-host-reagent.snakefile
ADD 06.map-reference-genome.snakefile /06.map-reference-genome.snakefile
ADD 06.bam_to_fastq.snakefile /06.bam_to_fastq.snakefile
ADD 07.map-by-gene.snakefile /07.map-by-gene.snakefile
ADD 08.call-variants-by-gene.snakefile /08.call-variants-by-gene.snakefile
ADD 09.minhash.dataset.snakefile /07.minhash.dataset.snakefile
ADD 09.minhash.sequences.snakefile /07.minhash.sequences.snakefile
ADD 22592-reagent-db.fasta.gz /22592-reagent-db.fasta.gz
ADD nt/ /nt 
