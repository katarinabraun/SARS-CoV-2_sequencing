#!/usr/bin/env python3
#Sniffles2
#Authors: Joseph Lalli and Kelsey Florek
#Performs SNP analysis of influenza genomes from raw reads.
import yaml
import argparse
import os, sys
import time, datetime
from shutil import copyfile
import multiprocessing as mp
from time import strftime
import glob
from operator import itemgetter

sys.path.append(os.path.join(os.getcwd(), 'supporting_code'))

import sc
from trim import trim
from mapping import mapping,indexing,average_depth
from consensus import consensus
import readcleaning as rc
from snpcaller import snpcaller
from RePlow import RePlow
from VCF_Results_Compiler import compareVCFs
from VCFannotater import VCFannotator
from SNPGenier import SNPgenier
import fileparser as fp
from send_notification_slack import sendSlack


#print main display title
sc.mainTitle()

os.environ["DOCKER_CLIENT_TIMEOUT"] = "240"
os.environ["COMPOSE_HTTP_TIMEOUT"]="240"

#determine command line arguments and get path
parser = argparse.ArgumentParser(description='Pipeline to examine SNPs from raw illumina reads')
parser.add_argument('-c',metavar='config', type=str,help="Location of configuration file",default='config.yml')
parser.add_argument('-i',metavar='input', type=str, help="Input fasta files directory - defaults to working directory")
parser.add_argument('-o',metavar='output', type=str,help="Output directory - defaults to sniffles_files")
parser.add_argument('-t',metavar='threads', type=int,help="number of cpus to use for pipeline",default=mp.cpu_count())
parser.add_argument('--debug',metavar='debug mode', type=str,help='Step to start analysis on. Options: trim, initial_mapping, indexing, norm_coverage, consensus, snpcalling, annotation, snpgenie')
parser.add_argument('--slack', metavar='slack user', type=str, help='Username of slack account to update during long runs. If left blank, won\'t send any slack messages.', default = False)
parser.add_argument('-B', action='store_true')#, help='The -B flag tells sniffles that it should process all bamfiles in the input directory, and jump directly to snpcalling')

if len(sys.argv[1:]) == 0:
		parser.print_help()
		parser.exit()
args = parser.parse_args()
numThreads = args.t
configFile = args.c
if args.debug:
	startingstep = ["trim", "initial_mapping", "indexing", "norm_coverage", "consensus", "snpcalling", "annotation", "snpgenie"].index(args.debug)
else:
	startingstep = 0

#trim, initial_mapping, indexing, annotation, snpgenier not supported as starting step

#get start time
start = time.time()


#open config file and store configuation
with open(configFile,'r') as ymlFile:
	cfg = yaml.safe_load(ymlFile)

cfg['slackUser'] = args.slack

#get input path
try:
	inDir = os.path.abspath(args.i)
except (AttributeError, TypeError) as err:
	inDir = os.getcwd()
	print(f"Raw reads directory {args.i} cannot be found. Sniffles will look for fasta files in the current working directory {inDir}.\n")

#create outdir
try:
	outDir = os.path.abspath(args.o)
except (AttributeError, TypeError) as err:
	outDir = os.getcwd()
	print(f"Output directory {args.o} cannot be found. Output will be placed in a separate folder in the current working directory {outDir}.\n")

sc.checkexists(outDir)

cfg['exec']['outdir'] = os.path.join(outDir,cfg['exec']['outdir'])

try:
	os.mkdir(cfg['exec']['outdir'])
except FileExistsError:
	cfg['exec']['outdir'] = cfg['exec']['outdir']+'_'+str(int(time.time()))
	os.mkdir(cfg['exec']['outdir'])
outDir = cfg['exec']['outdir']

logfile=os.path.join(outDir,cfg['exec']['logfile'])
cfg['exec']['logfile'] = logfile

startRunMessage = f"Beginning run at {strftime('%a, %d %b %Y %I:%M:%S %p', time.localtime())}"
sc.procTitle(startRunMessage, cfg)

with open(logfile,'a') as outlog:
	outlog.write(startRunMessage + "\n")

cfg['Errors'] = []

for reference, gtf in zip(list(cfg['exec']['referenceSequences']), cfg['postprocessing']['gtfFileNames']):
	sc.procTitle (f"Processing samples for reference sequence {reference.split('.')[0]}", cfg)
		
	#assign reference sequence value and gtf value to current reference and gtf
	cfg['exec']['referenceSequence'] = os.path.join(os.getcwd(), reference)
	cfg['postprocessing']['gtfFileName'] = os.path.join(os.getcwd(), gtf)

	try:
		inDir = os.path.abspath(args.i)
	except (AttributeError, TypeError) as err:
		inDir = os.getcwd()

	#check that fasta files are arranged appropriately
	if type(cfg['exec']['referenceSequences']) == list	:
		inDir = os.path.join(inDir, reference.split(".")[0])

		assert os.path.exists(inDir), "Error: Fasta files must be in subfolders named with the sequence to be referenced against"

		#create reference-specific outdir
		cfg['exec']['outdir'] = os.path.join(outDir, reference.split(".")[0])
		try:
			os.mkdir(cfg['exec']['outdir'])
		except FileExistsError:
			cfg['exec']['outdir'] = cfg['exec']['outdir']+'_'+str(int(time.time()))
			os.mkdir(cfg['exec']['outdir'])
	if not args.B:
		#copy reference to outdir
		copyfile(os.path.join(os.getcwd(), cfg['exec']['referenceSequence']),os.path.join(cfg['exec']['outdir'],reference))

		#trim the reads
		if startingstep <= 0:
			#parse and store read information from input directory
			readData = fp.RunFiles(inDir)
			ids = readData.return_id_list()
			trim(readData,cfg,numThreads)
			assert len(readData.runtime['trimmed']) > 0, "Cannot map reads: read trimming failed"
			#setup initial mapping jobs
			mapping_list = []
			print (readData.runtime['trimmed'])
			for id in readData.runtime['trimmed']:
				mapping_list.append((id,readData.runtime['trimmed'][id][0],readData.runtime['trimmed'][id][1],os.path.abspath(cfg['exec']['referenceSequence'])))
			
			if len(mapping_list) < len (ids):
				idsDidNotMakeIt = [id for id in ids if id not in list(map(itemgetter(1).replace('.fastq.gz',''), mapping_list))]
				cfg['Errors'].append(f'{", ".join(idsDidNotMakeIt)} were not successfully trimmed during run.')
		else:
			mapping_list = glob.glob(inDir+"trimmed/*.fastq.gz")
		
		
		sc.checkexists(os.path.join(cfg['exec']['outdir']+'/initial_mapping'))

		#index reference sequence
		#if startingstep <= 1:
		indexing(cfg,os.path.abspath(cfg['exec']['referenceSequence']))
		
		#run initial mapping jobs
		if startingstep <= 2:
			sc.procTitle("Beginning initial mapping", cfg)
			bam_list = mapping(cfg,mapping_list,cfg['exec']['outdir']+'/initial_mapping',numThreads)
			if len(bam_list) < len (mapping_list):
				idsDidNotMakeIt = [id for id in mapping_list if id not in [mapped.replace('.bam','.fastq.gz') for mapped in mapping_list]]
				cfg['Errors'].append(f'{", ".join(idsDidNotMakeIt)} were not successfully mapped during first mapping of run.')
		else:
			bam_list = glob.glob(inDir+"/initial_mapping/*.bam")

		

		#determine average depth
		

		if startingstep <= 2:
			#Filter by average depth coverage
			if cfg['exec']['coverageFilter']:
				origbamlist = bam_list
				assert len(bam_list) >0, "Cannot normalize coverage: 'mapping' failed to produce a list of bam files."
				bam_list = average_depth(cfg,bam_list,cfg['exec']['outdir']+'/initial_mapping',cfg['exec']['outdir']+'/coverage')
				if len(bam_list) < len (origbamlist):
					idsDidNotMakeIt = [id for id in origbamlist if id not in bam_list]
					cfg['Errors'].append(f'{", ".join(idsDidNotMakeIt)} did not meet coverage requirements.')

		#normalize coverage
		if startingstep <=3:
			if cfg['exec']['normalizeCoverage']:
				assert len(bam_list) >0, sc.procTitle("Cannot normalize coverage: 'average_depth' failed to produce a list of bam files.",cfg)
				fastq_list = rc.normCoverage(cfg,bam_list,numThreads)

				if not cfg['exec']['mapToConsensus']:
					mapping_list = []
					for fastq in fastq_list:
						read1 = fastq[0]
						read2 = fastq[1]
						if cfg['exec']['unpaired']:
							unpairedReads = fastq[2]
						else:
							unpairedReads = ""
						id = "_".join(os.path.basename(read1).split('_')[0:-1])
						mapping_list.append((id,read1,read2,os.path.abspath(cfg['exec']['referenceSequence']),unpairedReads))
					sc.checkexists(os.path.join(cfg['exec']['outdir']+'/norm_mapping'))
					origbamlist = bam_list
					bam_list = mapping(cfg,mapping_list,cfg['exec']['outdir']+'/norm_mapping',numThreads)
					if len(bam_list) < len (origbamlist):
						idsDidNotMakeIt = [id for id in origbamlist if id not in bam_list]
						cfg['Errors'].append(f'{", ".join(idsDidNotMakeIt)} were not successfully downsampled.')

		else:
			if cfg['exec']['mapToConsensus']:
				fwdfastqs = glob.glob(inDir+"/normalized/*_1.fastq.gz")
				revfastqs = glob.glob(inDir+"/normalized/*_2.fastq.gz")
				if cfg['exec']['unpaired']:
					unpairfastqs = glob.glob(inDir+"/normalized/*_U.fastq.gz")
					fastq_list = zip(fwdfastqs, revfastqs, unpairfastqs)
				else:
					fastq_list = zip(fwdfastqs, revfastqs)
			else:
				bam_list = glob.glob(inDir+"/norm_mapping/*.bam")


		if startingstep <= 4:
			#generate consensus
			if cfg['exec']['generateConsensus']:
				consensusRef = consensus(cfg,bam_list,numThreads)
				fasta_list = list(consensusRef)
				#fasta_list.append(consensusRef)
				# map reads to consensus
				if cfg['exec']['generateConsensus'] and cfg['exec']['mapToConsensus']:
					mapping_list = []
					indexing(cfg,*fasta_list)
					for fastq in fastq_list:
						print (fastq)
						read1 = fastq[0]
						read2 = fastq[1]
						if cfg['exec']['unpaired']:
							unpairedReads = fastq[2]
						else:
							unpairedReads = ""
						id = "_".join(os.path.basename(read1).split('_')[0:-1])
						mapping_list.append((id,read1,read2,consensusRef,unpairedReads))
					
					sc.checkexists(os.path.join(cfg['exec']['outdir']+'/map_to_consensus'))
					origbamlist = bam_list
					bam_list = mapping(cfg,mapping_list,cfg['exec']['outdir']+'/map_to_consensus',numThreads)
					if len(bam_list) < len (origbamlist):
						idsDidNotMakeIt = [id for id in origbamlist if id not in bam_list]
						cfg['Errors'].append(f'{", ".join(idsDidNotMakeIt)} were not successfully mapped during the second mapping of the run.')
					cfg['exec']['referenceSequence'] = consensusRef
				elif cfg['exec']['generateConsensus']:
					print ("You have told Sniffles to map your sequences to your sample reference, but have not selected the \"Generate Consensus\" option in your config file. There is no sample consensus to map against. Sniffles will map samples to the reference provided.")
		elif cfg['exec']['mapToConsensus']:
			bam_list = glob.glob(inDir+"/map_to_consensus/*.bam")
			cfg['exec']['referenceSequence'] = os.path.join(inDir, "ref_sequence", os.path.basename(cfg['exec']['referenceSequence']).split(".")[0]+"_consensus_noambig.fasta")

		#call snps
		print (bam_list)
		assert len(bam_list) > 0, "Cannot call SNPs: Mapping reads to consensus failed."
	else: #(if args.B:)
		indexing(cfg,os.path.abspath(cfg['exec']['referenceSequence']))
		bam_list = glob.glob(inDir+'/*.bam')

	if startingstep <= 5:
		if cfg['exec']['replicates']:
			if cfg['exec']['callSNPs'] == "Varscan":
				snpcaller(cfg,bam_list,numThreads)
			elif cfg['exec']['callSNPs'] == "RePlow":
				RePlow(cfg,bam_list,numThreads)
			elif cfg['exec']['callSNPs'] == "Compare":
				replow = RePlow(cfg,bam_list,numThreads)
				varscan = snpcaller(cfg,bam_list,numThreads)
				compareVCFs(varscan, replow) #creates merged vcf that compares the output of the two SNP callers
			else:
				raise ValueError('callSNPs in config.yml must be either \'Varscan\' or \'RePlow\'')
		else:
			snpcaller(cfg,bam_list,numThreads)
		vcfpath = os.path.join(cfg['exec']['outdir'],"snp_calls")
		vcfs = glob.glob(vcfpath + "/[0-9].vcf")
		if len(vcfs) < len(bam_list):
			idsDidNotMakeIt = [id.split('/')[-1] for id in bam_list if id not in [vcf.replace('.vcf','.bam') for vcf in vcfs]]
			cfg['Errors'].append(f'{", ".join(idsDidNotMakeIt)} was/were not successfully snp called during run.')

	else:
		vcfpath = os.path.join(inDir,"snp_calls")
		print (f"vcfpath: {vcfpath}")
	
	if startingstep <= 6:
		if cfg['exec']['annotateSNPs']:
			VCFannotator(cfg, vcfpath)

	if startingstep <=7:
		if cfg['exec']['SNPgenier']:
			SNPgenier(cfg, vcfpath, numThreads)
	

end = time.time()
runtime = round(end - start,2)
runtime = str(datetime.timedelta(seconds=runtime))
sc.procTitle(f'Sniffles: Finished with a total runtime of {runtime}.', cfg)

for error in cfg['Errors']:
	print (error)

with open(logfile,'a') as outlog:
	outlog.write('\n' + '\n'.join(cfg['Errors']))
