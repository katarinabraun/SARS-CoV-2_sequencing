# parameters
CLEANED_FASTQ_PATH = config['cleaned']
NT_SKETCH_PATH = '/nt/global/projectb/sandbox/gaag/bbtools/nt/current/*.sketch'

# get list of every FASTQ file
ONT_FASTQ, = glob_wildcards(CLEANED_FASTQ_PATH  + '/{barcode}.clean.fastq.gz')

rule all:
    input:
        expand('minhash/{id}.sequences.txt', id = ONT_FASTQ),

rule minhash_sequences:
    input:
        CLEANED_FASTQ_PATH + '/{id}.clean.fastq.gz'
    output:
        'minhash/{id}.sequences.txt'
    run:
        shell('comparesketch.sh \
            qin=33 \
            in={input[0]} \
            out={output[0]}\
            mode=sequence \
            size=100000 \
            minhits=1 \
            ow=t \
            ' + NT_SKETCH_PATH)

    