# parameters
MAPPED_PATH = config['mapped']
MAPPING_GENOME = config['mapping_genome']

# get list of every FASTQ file
ONT_FASTQ, = glob_wildcards(MAPPED_PATH  + '/{barcode}.primary.bam')

rule all:
    input:
        expand('mapped/{id}.fastq', id = ONT_FASTQ),

rule convert_bam_to_fastq:
    input:
        MAPPED_PATH + '/{id}.primary.bam',
        MAPPING_GENOME
    output:
        'mapped/{id}.fastq'
    run:
        shell('reformat.sh \
            in={input[0]} \
            out={output[0]} \
            ref={input[1]} \
            mappedonly \
            ow=t')
