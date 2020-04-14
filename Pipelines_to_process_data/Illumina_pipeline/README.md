# Sniffles -- two versions 

There are two versions of the pipeline. 

`from_bams` takes in bam files (alignment files) as input -- we recommend using this pipeline to replicate our analyses. The input files for this pipeline can be found in `data_raw/Illumina_bams/*`. 

`from_fastQs` takes in paired fastQ files (R1 and R2). We have included this version of the pipeline for thoroughness, but we do not provide the raw, unpaired, unmerged Illumina fastQ files on the SRA bioproject or within this repository because we wanted to ensure we were not uploading any contaminating human reads. 