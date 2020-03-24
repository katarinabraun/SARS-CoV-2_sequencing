## parameters ##
import ntpath 

# expects demultiplexed FASTQ
BYGENEBAM = config['bygenebam']
ONT_BASENAME = ntpath.basename(BYGENEBAM).replace('.primary.bam', '')

# reference databases
MAPPING_GENOME = config['mapping_genome']

CALLVARIANTS = 'callvariants.sh \
in={input[0]} \
ref={input[1]} \
minallelefraction=0.10 \
rarity=0.10 \
calldel=f \
callindel=f \
out={output[0]}'

## rules
rule all:
    input:
        expand('mapped_bygene/{id}.vcf', id = ONT_BASENAME)
        
rule callvariants:
	input:
		BYGENEBAM,
		MAPPING_GENOME,
	output:
		'mapped_bygene/{id}.vcf'
	run:
		shell(CALLVARIANTS)