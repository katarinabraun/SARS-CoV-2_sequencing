# parameters
CLEANED_FASTQ_PATH = config['merged_demultiplexed']

# get list of every FASTQ file
ONT_FASTQ, = glob_wildcards(CLEANED_FASTQ_PATH  + '/{barcode}.fastq.gz')

rule all:
    input:
        expand('subsample/{id}.fastq.gz', id = ONT_FASTQ),

rule subsample_QC:
    input:
        CLEANED_FASTQ_PATH + '/{id}.fastq.gz'
    output:
        'subsample/{id}.fastq.gz'
    run:
        shell('reformat.sh \
            in={input[0]} \
            out={output[0]} \
            forcetrimleft=30 \
            forcetrimright2=30 \
            mincalledquality=7 \
            ow=t')