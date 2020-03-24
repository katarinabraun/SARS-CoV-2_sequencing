# parameters
CLEANED_FASTQ_PATH = config['cleaned']
NT_SKETCH_PATH = '/nt/global/projectb/sandbox/gaag/bbtools/nt/current/*.sketch'

# get list of every FASTQ file
ONT_FASTQ, = glob_wildcards(CLEANED_FASTQ_PATH  + '/{barcode}.clean.fastq.gz')

rule all:
    input:
        expand('minhash/{id}.dataset.txt', id = ONT_FASTQ),

rule minhash_dataset:
    input:
        CLEANED_FASTQ_PATH + '/{id}.clean.fastq.gz'
    output:
        'minhash/{id}.dataset.txt'
    run:
        shell('comparesketch.sh \
            qin=33 \
            in={input[0]} \
            out={output[0]}\
            size=100000 \
            records=1000 \
            minhits=1 \
            ow=t \
            ' + NT_SKETCH_PATH)
