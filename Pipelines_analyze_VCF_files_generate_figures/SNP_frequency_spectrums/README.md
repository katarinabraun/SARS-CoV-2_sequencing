# SNV Frequency Spectrum (SFS)

**Author**: Katarina Braun 

## Overview
This notebook plots SNV frequency spectrums. It pulls SNVs and their frequencies from the `-cleaned.csv` files that were generated in the `SNVs.ipynb` script. 

## Dependencies 
- pandas
- numpy
- matplotlib
- glob
- os
- functools
- pandas
- itertools
- glob
- os
- sklearn 

### A useful description regarding the interpretation of these analyses 
**From Moncla et al, 2019, bioRxiv (https://doi.org/10.1101/683151):**  
"Purifying selection removes new variants from the population, generating an excess of low-freq variants, while positive selection promotes accumulation of high-frequency polymorphisms. Exponential population expansion also causes excess low-frequency variation; however, while selection disproportionately affects nonsynonymous variants, demographic factors affect synonymous and nonsynonymous variants equally."

**The figure will look something like this**: 
- Y = proportion of SNVs
- X = within-host SNV frequency bins: 1-10%, 10-20%, 20-30%, 30-40%, 40-50%

I will also try to derive the "neutral expectation" -- that is the distribution of SNPs expected for a given population assuming that the population is not under selection and is at some sort of steady-state equilibrium. 

For the neutral expectation, Trevor Bedford suggests this will follow a 1/x distribution. I can then just integrate over a 1/x distribution between each bin size (0.01 to 0.1, 0.1 to 0.2, etc...). Then I'll calculate the proportion of the total that fall into each bin. Dr. Louise Moncla (cited above) already prepared a notebook to do this: `neutral-expectation.ipynb`, which I am going to utilize here. 
    
## Input: 

`VCFs_Illumina/cleaned/*`

## Output: 

`SFS.pdf`