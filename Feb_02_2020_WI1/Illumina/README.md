# Sniffles

Sniffles is a pipeline for the analysis of influenza genomes. It's a configurable pipeline that performs multiple functions to generate variant and consensus level information. It requires minimal dependencies as the pipeline relies on a docker container to host the software.

## Table of Contents
* [Requirements](#requirements)
* [Installing](#installing)
* [Usage](#usage)

### Requirements
* Linux or MacOS
* Python 3.6 or later
* Docker

### Installing
* Make sure you are using the correct version of python
* Use git to get the latest sniffles build `git clone https://github.com/k-florek/sniffles.git`
* Install the [Docker CE engine](https://docs.docker.com/install/)
* Install the required python libraries by pointing the pip installer to the requirements document in the sniffles project `pip3 install -r requirements.txt`

### Usage
Sniffles uses a configuration file to provide parameters to the pipeline. There is a default configuration file included in the sniffles project. To run the pipeline simple point the `sniffles.py` program at the correct config.yml and directory containing the raw Illumina reads.

```
 #####
#     #  #    #  #  ######  ######  #       ######   ####
#        ##   #  #  #       #       #       #       #
 #####   # #  #  #  #####   #####   #       #####    ####
      #  #  # #  #  #       #       #       #            #
#     #  #   ##  #  #       #       #       #       #    #
 #####   #    #  #  #       #       ######  ######   ####



usage: sniffles.py [-h] [-c config] [-i input] [-o output] [-t threads]

Pipeline to examine SNPs from raw illumina reads

optional arguments:
 -h, --help  show this help message and exit
 -c config   config file
 -i input    raw reads directory - defaults to working directory
 -o output   output directory - defaults to working directory
 -t threads  number of cpus to use for pipeline
```

Example usage:
`./sniffles.py -c config.yml -i ~/my_reads/ -t 8 -o my_results`
