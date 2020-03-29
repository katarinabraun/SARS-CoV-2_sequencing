### Authors

**Katarina Braun** (graduate student in Thomas Friedrich's lab)  Twitter: [@KATarinambraun](https://twitter.com/KATarinambraun) GitHub: [katarinabraun](https://github.com/katarinabraun)

**Gage Moreno** (graduate student in Dave O'Connor's lab) Twitter: [@GageKMoreno](https://twitter.com/GageKMoreno)

---------------------
## Data availability

We have made our cleaned FASTQ* files available at the following SRA: [PRJNA614504]()

**Note these FASTQs have been depleted of host sequences and other contaminating sequences.

These analyses and scripts are additionally available on GitHub [SARS-CoV-2_sequencing](https://github.com/katarinabraun/SARS-CoV-2_sequencing)

---------------------
## Scripts used
### ONT scripts  
These data were analyzed to produce consensus sequences and to generate VCFs with consensus-level SNVs using the [ARTIC bioinformatic pipeline](https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html). 

We used a custom script to do the ARTIC bioinformatic pipeline in real time after a barcoded sample reaches 100k reads. This will allow the user to only run a flow cell for a few hours and collect consensus data in real time. This script is uploaded as `artic_rt_workflow.sh` but note that it needs to be configured to work with your desired network. Briefly, this script will watch the fastq_pass folder being made on the GridION and demultiplex each fastq file generated in real time. It will then count the number of barcodes found. Once a barcode reaches 100k reads, it will trigger the rest of the ARTIC bioinformatics workflow which will map to the severe acute respiratory syndrome coronavirus isolation from Wuhan, China (Genbank: MN908947.3) using `mimimap2`. This alignment will then be used to generate consensus sequences and variant calls.
- Note: Configure the bash script at lines 28 and 42. 

Samples were visualized using [RAMPART](https://artic-network.github.io/rampart/)

## Methods
We sequenced 15 positive samples from Madison, WI. We started with 3µl of nasal swab sample diluted in 197µl of water. This 200µl original dilution was then processed to extract viral RNA. Two of these samples had very viral loads by QRT-PCR and are still on the GridION because we have not achieved consensus-level data yet. We will update this page when we have full coverage available on those two additional samples. 

### Sample information 

| Sample                 | GISAID                                                        | SRA Accession | SRR Run # | BioSample | BioProject  | Collection Date |
|------------------------|---------------------------------------------------------------|---------------|-----------|-----------|-------------|-----------------|
| hCoV-19/USA/WI-11/2020 | [EPI_ISL_417505](https://www.epicov.org/epi3/frontend#516793) |               |           |           | PRJNA614504 | Unknown         |
| hCoV-19/USA/WI-12/2020 | [EPI_ISL_417506](https://www.epicov.org/epi3/frontend#3ff138) |               |           |           | PRJNA614504 | Unknown         |
| hCoV-19/USA/WI-13/2020 | [EPI_ISL_417516](https://www.epicov.org/epi3/frontend#3c5e53) |               |           |           | PRJNA614504 | Unknown         |
| hCoV-19/USA/WI-14/2020 | [EPI_ISL_417513](https://www.epicov.org/epi3/frontend#1eeaad) |               |           |           | PRJNA614504 | Unknown         |
| hCoV-19/USA/WI-15/2020 | [EPI_ISL_417504](https://www.epicov.org/epi3/frontend#52f026) |               |           |           | PRJNA614504 | Unknown         |
| hCoV-19/USA/WI-16/2020 | [EPI_ISL_417509](https://www.epicov.org/epi3/frontend#4233a4) |               |           |           | PRJNA614504 | Unknown         |
| hCoV-19/USA/WI-17/2020 | [EPI_ISL_417517](https://www.epicov.org/epi3/frontend#514901) |               |           |           | PRJNA614504 | Unknown         |
| hCoV-19/USA/WI-18/2020 | [EPI_ISL_417515](https://www.epicov.org/epi3/frontend#4fc3b9) |               |           |           | PRJNA614504 | Unknown         |
| hCoV-19/USA/WI-19/2020 | [EPI_ISL_417510](https://www.epicov.org/epi3/frontend#d6c50)  |               |           |           | PRJNA614504 | Unknown         |
| hCoV-19/USA/WI-20/2020 |                                                               |               |           |           |             | Unknown         |
| hCoV-19/USA/WI-21/2020 | [EPI_ISL_417508](https://www.epicov.org/epi3/frontend#5d9370) |               |           |           | PRJNA614504 | Unknown         |
| hCoV-19/USA/WI-22/2020 | [EPI_ISL_417514](https://www.epicov.org/epi3/frontend#622fb6) |               |           |           | PRJNA614504 | Unknown         |
| hCoV-19/USA/WI-23/2020 | [EPI_ISL_417507](https://www.epicov.org/epi3/frontend#4eebd2) |               |           |           | PRJNA614504 | Unknown         |
| hCoV-19/USA/WI-24/2020 | [EPI_ISL_417512](https://www.epicov.org/epi3/frontend#4dea1d) |               |           |           | PRJNA614504 | Unknown         |

** hCoV-19/USA/WI-20/2020 will be added at a later time 

### Sequencing and analysis methods 
We isolated viral RNA using the Maxwell 48. Complete protocol can be found [here](https://openresearch.labkey.com/wiki/Coven/download.view?entityId=bf8f06ee-501b-1038-80cc-8e513fde0084&name=Maxwell%20RSC%20Viral%20Total%20Nucleic%20Acid%20Purification%20Kit%20TM420.pdf)

We prepared each sample for ONT sequencing using the ARTIC network tiled amplicon approach: 
ARTIC [complete protocol](https://www.protocols.io/view/ncov-2019-sequencing-protocol-bbmuik6w).

We then prepared libraries for Oxford Nanopore Sequencing on a GridION using the 1D ligation sequencing kit (SQK-LSK109) with its native barcodes and sequenced on an R9.4 flow cell. 

All reads were processed using the `artic_rt_workflow.sh` script. Briefly, this script will watch the fastq_pass folder being made on the GridION and demultiplex each fastq file generated in real time. It will then count the number of barcodes found. Once a barcode reaches 100k reads, it will trigger the rest of the ARTIC bioinformatics workflow which will map to the severe acute respiratory syndrome coronavirus isolation from Wuhan, China (Genbank: MN908947.3) using `mimimap2`. This alignment will then be used to generate consensus sequences and variant calls. 

## Results

#### Coverage across SARS-CoV-2 genome  
We pulled data off of the GridION after we generated 100K reads for each sample, averaging to 400x coverage across the entire genome. We are not sequencing beyond this depth right now because we are most interested in reporting/uploading consensus-level sequences to GISAID and Nextstrain. 


#### Consensus level differences between samples

We mapped all swabs against the Wuhan reference sequence (MN908947.3) and identified consensus-level differences from the reference sequence and the thirteen swabs we sequenced here. We also took a quick look at the impact of each of these SNVs and have included these annotations in the second table. 

![Here’s a table outlining consensus-level SNV differences among samples.](consensus_SNV_diferences.png)

Overall, consensus-level SNVs are most likely to fall in ORF1ab (which occupies the majority of the genome). We have identified very few SNVs in the S gene, which encodes Spike protein expressed on the surface of the virion -- and is therefore antigenically most relevant. 
There are two few consensus-level SNVs to make confident inferences regarding population selective pressures. We have not called minor variants at this time, but hope to do this soon. 

#### Nextstrain phylogenies 
Consensus sequences were uploaded to GISAID and subsequently added to Nextstrain for phylogenetic analyses. Sample naming scheme is in the table below. 

| Sample                 | Nextstain                                                              | GISAID                                                        |
|------------------------|------------------------------------------------------------------------|---------------------------------------------------------------|
| hCoV-19/USA/WI-11/2020 | [hCoV-19/USA/WI-11/2020](https://nextstrain.org/ncov?s=USA/WI-11/2020) | [EPI_ISL_417505](https://www.epicov.org/epi3/frontend#516793) |
| hCoV-19/USA/WI-12/2020 | [hCoV-19/USA/WI-12/2020](https://nextstrain.org/ncov?s=USA/WI-12/2020) | [EPI_ISL_417506](https://www.epicov.org/epi3/frontend#3ff138) |
| hCoV-19/USA/WI-13/2020 | [hCoV-19/USA/WI-13/2020](https://nextstrain.org/ncov?s=USA/WI-13/2020) | [EPI_ISL_417516](https://www.epicov.org/epi3/frontend#3c5e53) |
| hCoV-19/USA/WI-14/2020 | [hCoV-19/USA/WI-14/2020](https://nextstrain.org/ncov?s=USA/WI-14/2020) | [EPI_ISL_417513](https://www.epicov.org/epi3/frontend#1eeaad) |
| hCoV-19/USA/WI-15/2020 | [hCoV-19/USA/WI-15/2020](https://nextstrain.org/ncov?s=USA/WI-15/2020) | [EPI_ISL_417504](https://www.epicov.org/epi3/frontend#52f026) |
| hCoV-19/USA/WI-16/2020 | [hCoV-19/USA/WI-16/2020](https://nextstrain.org/ncov?s=USA/WI-16/2020) | [EPI_ISL_417509](https://www.epicov.org/epi3/frontend#4233a4) |
| hCoV-19/USA/WI-17/2020 | [hCoV-19/USA/WI-17/2020](https://nextstrain.org/ncov?s=USA/WI-17/2020) | [EPI_ISL_417517](https://www.epicov.org/epi3/frontend#514901) |
| hCoV-19/USA/WI-18/2020 | [hCoV-19/USA/WI-18/2020](https://nextstrain.org/ncov?s=USA/WI-18/2020) | [EPI_ISL_417515](https://www.epicov.org/epi3/frontend#4fc3b9) |
| hCoV-19/USA/WI-19/2020 | [hCoV-19/USA/WI-19/2020](https://nextstrain.org/ncov?s=USA/WI-19/2020) | [EPI_ISL_417510](https://www.epicov.org/epi3/frontend#d6c50)  |
| hCoV-19/USA/WI-20/2020 | [hCoV-19/USA/WI-20/2020](https://nextstrain.org/ncov?s=USA/WI-20/2020) |                                                               |
| hCoV-19/USA/WI-21/2020 | [hCoV-19/USA/WI-21/2020](https://nextstrain.org/ncov?s=USA/WI-21/2020) | [EPI_ISL_417508](https://www.epicov.org/epi3/frontend#5d9370) |
| hCoV-19/USA/WI-22/2020 | [hCoV-19/USA/WI-22/2020](https://nextstrain.org/ncov?s=USA/WI-22/2020) | [EPI_ISL_417514](https://www.epicov.org/epi3/frontend#622fb6) |
| hCoV-19/USA/WI-23/2020 | [hCoV-19/USA/WI-23/2020](https://nextstrain.org/ncov?s=USA/WI-23/2020) | [EPI_ISL_417507](https://www.epicov.org/epi3/frontend#4eebd2) |
| hCoV-19/USA/WI-24/2020 | [hCoV-19/USA/WI-24/2020](https://nextstrain.org/ncov?s=USA/WI-24/2020) | [EPI_ISL_417512](https://www.epicov.org/epi3/frontend#4dea1d) |

#### Preliminary analysis of transmission clusters 

| Sample Name            | Clade Ident. | Proximal geographic origin | Transmission Chain  |
|------------------------|--------------|----------------------------|---------------------|
| hCoV-19/USA/WI-11/2020 | A2a          | Europe                     | WI4                 |
| hCoV-19/USA/WI-12/2020 | A2a          | Europe                     | **distinct (WI12)**     |
| hCoV-19/USA/WI-13/2020 | A2a          | Europe                     | **distinct (WI13)**     |
| hCoV-19/USA/WI-14/2020 | A2a          | Europe                     | WI12                |
| hCoV-19/USA/WI-15/2020 | B1           | WA                         | WI9                 |
| hCoV-19/USA/WI-16/2020 | A2a          | Europe                     | **distinct (WI16)**     |
| hCoV-19/USA/WI-17/2020 | A2a          | Europe via DR Congo        | **distinct (WI17)**     |
| hCoV-19/USA/WI-18/2020 | A2a          | Europe                     | **distinct (WI8)**      |
| hCoV-19/USA/WI-19/2020 | A2a          | Europe                     | WI12                |
| hCoV-19/USA/WI-20/2020 |              |                            |                     |
| hCoV-19/USA/WI-21/2020 | A2a          | Europe                     | WI12                |
| hCoV-19/USA/WI-22/2020 | B4           | Australia                  | WI6                 |
| hCoV-19/USA/WI-23/2020 | B1           | WA                         | **distinct (WI23)**     |
| hCoV-19/USA/WI-24/2020 | A2a          | Europe                     | **distinct (WI16)**     |


Based on Nextstrain phylogenetic analysis, it looks like there may be ≥ 10 separate introductions of SARS-CoV-2 into Wisconsin. This is very preliminary analysis and we are working on more sophisticated analysis regarding the number of unique introductions into Wisconsin now.