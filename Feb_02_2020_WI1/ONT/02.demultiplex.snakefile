# parameters
MIN_LENGTH = config['min_length']
PARTITIONED_FASTQ_PATH = config['partitioned_fastq_folder']

# get list of every FASTQ file
ONT_FASTQ, = glob_wildcards(PARTITIONED_FASTQ_PATH  + '/{partition}.fastq')

rule all:
    input:
        # create input file after running qcat on each sample
        expand(PARTITIONED_FASTQ_PATH + '/{id}.txt', id = ONT_FASTQ)

rule qcat:
    '''run qcat to trim barcodes and adapters from ONT fastq files
    use a single thread for each qcat operation so each file gets processed separately'''
    input: 
        PARTITIONED_FASTQ_PATH + '/{id}.fastq'
    output:
        PARTITIONED_FASTQ_PATH + '/{id}.txt'
    params:
        ID = '{id}'
    threads: 1
    run:
        shell('mkdir -p demultiplexed/{params.ID} && \
            qcat \
            -f {input[0]} \
            -b demultiplexed/{params.ID} \
            --min-read-length ' + str(MIN_LENGTH) + ' \
            --trim && \
            touch {output}')