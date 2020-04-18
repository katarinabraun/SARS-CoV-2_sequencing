# ARTIC Analysis Pipeline 

We used a custom script to do the ARTIC bioinformatic pipeline in real time after a barcoded sample reaches 100k reads. This will allow the user to only run a flow cell for a few hours and collect consensus data in real time. This script is uploaded as artic_rt_workflow.sh but note that it needs to be configured to work with your desired network. Briefly, this script will watch the fastq_pass folder being made on the GridION and demultiplex each fastq file generated in real time. It will then count the number of barcodes found. Once a barcode reaches 100k reads, it will trigger the rest of the ARTIC bioinformatics workflow which will map to the severe acute respiratory syndrome coronavirus isolation from Wuhan, China (Genbank: MN908947.3) using mimimap2. This alignment will then be used to generate consensus sequences and variant calls.

Note: Configure the bash script at lines 28 and 42.

Samples were visualized using [RAMPART](https://artic-network.github.io/rampart/).