#!/usr/bin/env python
# coding: utf-8

###Written by Joseph Lalli and Louise Moncla, PhD, edited by Kat Braun

import sys, glob, os, re
import pandas as pd
import numpy as np
from sc import procTitle,checkexists
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
import time, datetime
from tqdm import tqdm



# dictionary to convert single letter to 3 letter amino acid symbols
amino_acid_abbreviations = {"A": "Ala","R": "Arg","N": "Asn","D": "Asp","C": "Cys","Q": "Gln","E": "Glu","G": "Gly","H": "His","I": "Ile","L": "Leu","K": "Lys","M": "Met","F": "Phe","P": "Pro","O": "Pyl","S": "Ser","U": "Sec","T": "Thr","W": "Trp","Y": "Tyr","V": "Val", "B":"Asx","Z":"Glx","U":"Sec","X":"Xaa","J":"Xle","*":"Stop"}
start_codon = "atg"
stop_codons = ["taa","tag","tga"]
################################


def VCFannotator(runCFG, vcffiles):
	# read in reference sequences and gtfs and store information about protein sequences
	# create a dictionary of gene names and start/stop sites, allowing for more than one start/stop site.

	#import file location parameters from config file
	outDir = os.path.join(runCFG['exec']['outdir'], 'vcf_annotations')
	checkexists(outDir)
	logfile = os.path.join(outDir,runCFG['exec']['logfile'])
	refseqfasta = runCFG['exec']['referenceSequence']
	refseqname = refseqfasta.split(".")[0]

	if runCFG['exec']['mapToConsensus']:
		refseqfasta = os.path.join(runCFG['exec']['outdir'], 'ref_sequence', refseqfasta)

	#get start time
	start1 = time.time()

	procTitle('Annotating SNPs', runCFG)
	
	#Extract coding sequence coordinates from gtf files:
	coding_regions = {} #will be dictionary of dictionaryies (format segment:gene:[[startExon1, stopExon1], [startExon2, stopExon2]])
	with open(runCFG['postprocessing']['gtfFileName'], "r") as gtf:
		for line in gtf:
			if line.strip("\n") != "":	  # ignore blank lines (otherwise throws an index error)
				line = line.replace("/", "_")
				lineitems = line.split("\t")
				segment_name = lineitems[0]
				annotation_type = lineitems[2]
				start = int(lineitems[3]) - 1  # adding the -1 here for 0 indexing
				stop = int(lineitems[4]) - 1	# adding the -1 here for 0 indexing
				gene_name = lineitems[8]
				gene_name = gene_name.split(";")[0]
				gene_name = gene_name.replace("gene_id ","")
				gene_name = gene_name.replace("\"","")

				if annotation_type.lower() == "cds":
					if segment_name not in coding_regions:
						coding_regions[segment_name] = {}
						coding_regions[segment_name][gene_name] = [[start, stop]]
					elif segment_name in coding_regions and gene_name not in coding_regions[segment_name]:
						coding_regions[segment_name][gene_name] = [[start, stop]]
					elif gene_name in coding_regions[segment_name]:
						coding_regions[segment_name][gene_name].append([start, stop])
		
	# pull in reference fasta file, separate gene segments into a dictionary
	ref_segments = {}

	for seq in SeqIO.parse(refseqfasta, "fasta"):
		refseqname = str(seq.id).replace("/", "_")
		sequence = str(seq.seq).lower()
		ref_segments[refseqname] = sequence
		
	# use gene coordinates to create coding sequences from reference sequences
	transcripts = {}

	#Reminder of current data structures:
	#coding_regions[segment][gene]:coordinates of genes
	#ref_segments[nameofsegment]:sequence
	
	for segment in coding_regions:
		for gene in coding_regions[segment]:
			transcripts[gene] = ""
			coordinates = coding_regions[segment][gene]  # define the coding regions for each gene
			for start, stop in coordinates:   # loop through start/stop sites in coding regions
				sequence_chunk = ref_segments[segment][start:stop+1]
				transcripts[gene] = transcripts[gene] + sequence_chunk	 # append each piece of the transcript together 
	
	# loop through each transcript to make sure that it begins with a start codon and ends with a stop codon
	#for t in transcripts: 
		#if transcripts[t][0:3] != start_codon:
			#print("WARNING! " + refseqname + " " + t + " does not contain a start codon! The first three nucleotides are " + transcripts[t][0:3])
		#if transcripts[t][-3:] not in stop_codons:
			#print("WARNING! " + refseqname + " " + t + " does not contain a stop codon! These are the last 3 nucleotides: " + transcripts[t][-3:])


	print (vcffiles)
	if os.path.isdir(vcffiles):
		vcffiles = glob.glob(vcffiles + "/*.vcf")
	elif type(vcffiles) == list:
		if vcffiles[0].split(".")[-1] == "vcf":
			pass
	else:
		print ("vcffiles has no vcf files!")

	listofmutstoExport=[]
	##Loop through each vcf file and annotate amino acid changes
	print (vcffiles)
	for i, vcfname in tqdm(enumerate(vcffiles)):
		with open (vcfname, "r") as TextVCF:
			for index, line in enumerate(TextVCF, 0):
				if "#CHROM" in line:
					rowstoskip = index
		
		#Reads the vcf file into a pandas DataFrame
		print (vcfname)
		try:
			vcfDF = pd.read_csv(vcfname, sep='\t', skiprows=rowstoskip)
		except OSError as inst:
			print ("\n" + vcfname + " did not open appropriately. Please check file.\n")	
		
		#fix the / in chrome bug
		vcfDF["#CHROM"] = vcfDF["#CHROM"].str.replace("/", "_")
		#extract frequencies for list of muts to export:
		#In order to make this easier, I'm going to assume each VCF has only one sample. This code DOES NOT WORK for more than one sample per VCF.
		#print(vcfDF.iloc[:,-1].str.split(":").str[6].str.rstrip('%').astype('float')/100)
		freqlocation = vcfDF.loc[0, "FORMAT"].split(":").index("FREQ")
		try:
			vcfDF["FREQ"] = vcfDF.iloc[:,-1].str.split(":").str[freqlocation].astype('float')
		except ValueError:
			vcfDF["FREQ"] =  vcfDF.iloc[:,-1].str.split(":").str[freqlocation].str.rstrip('%').astype('float')/100
		except:
			raise

		listofmuts=[]
		
		#loop through each line in vcfDF, extract chrom, pos, reference nucleotide, alternate nucleotide
		for chrom, pos, ref, alt, freq in zip(vcfDF['#CHROM'], vcfDF["POS"], vcfDF["REF"].str.lower(), vcfDF["ALT"].str.lower(), vcfDF["FREQ"]):
			pos-=1 #subtract one from position to convert from VCF's 1 indexing to python's 0 
			for gene in coding_regions[chrom].keys(): #loop through each gene potentially applicable to that position (i.e., all on chromosome)
				priorExonLength=0
				#print (gene)
				for start, stop in coding_regions[chrom][gene]:    #loop through each exon of gene
					#if pos in exon, calculate codon, reference aa, and variant aa
					#print (f"pos: {pos} start: {start} stop: {stop}")
					if pos in range(start,stop):
						#print ('is in range, annotating.')
						within_gene_position = pos - start + priorExonLength #within gene position is the position in this exon (pos-startOfExon), plus the length of any prior exons (exonstart)
						codon_pos = (within_gene_position % 3)
						alternatetranscript = transcripts[gene][:within_gene_position]+ alt + transcripts[gene][within_gene_position+1:]
						
						codon = transcripts[gene][(within_gene_position - codon_pos):(within_gene_position+(3-codon_pos))]
						variantcodon = alternatetranscript[(within_gene_position - codon_pos): (within_gene_position+(3-codon_pos))]
						ref_aa = Seq(codon).translate()
						variant_aa = Seq(variantcodon).translate()
						aa_num = str(int(within_gene_position/3)+1)
						#Catch errors in annotation calculations where the math results in an incorrect codon
						if codon[codon_pos] != ref:
							print ("Something's quite wrong here. The reference SNP is not what it should be.")
							print (f"\n\nchrom: {chrom}, gene: {gene}")
							print(f"\npos: {pos}  within_gene_position: {within_gene_position}\ncodon_pos: {codon_pos}  codon: {codon}  variantcodon: {variantcodon}\n\n")
							print (f"ref:  {ref} alt: {alt}  ref_aa: {ref_aa} \nvariant_aa: {variant_aa}\n aa_num: {aa_num}\n")
							print (transcripts[gene])
						ref_aa = Seq(codon).translate()
						variant_aa = Seq(variantcodon).translate()
						aa_num = str(int(within_gene_position/3)+1)
						if ref_aa != variant_aa:
							listofmuts.append([chrom, gene, pos+1, str(ref_aa + aa_num + variant_aa)])
							if freq > 0.01 and freq < 0.99:
								listofmutstoExport.append({"segment": chrom, 'gene': gene, 'position': pos+1, 'frequency': float(freq), 'AAchange': str(ref_aa + aa_num + variant_aa)})
						elif ref_aa == variant_aa:
							listofmuts.append([chrom, gene, pos+1, "."])
						break #if pos is in exon, stop looping though exons
					else:
						priorExonLength += (stop+1-start) #The next exon will begin after the length of this exon, i.e., after the stop point minus the start point
				#else statement only executed if for loop finishes without breaking (ie if pos is never within the gene being examined)
				else: 
					listofmuts.append([chrom, gene, pos+1, "not in ORF"])
					continue #continue onto next gene

		AAchange = pd.DataFrame(listofmuts, columns=['#CHROM', 'gene', 'POS', 'AAchange'])

		vcfDF = vcfDF.merge(AAchange, how='left', on=['#CHROM', 'POS'])

		vcfDF['gene'] = vcfDF['gene'].astype(str)
		vcfDF['gene'] = vcfDF['gene'].replace("NA", "NA gene")
		annotatedVCFname = os.path.basename(vcfname).split(".")[0] + ".annotated_vcf"
		outputfile = vcfDF.to_csv(os.path.join(outDir, annotatedVCFname), sep = "\t", index = None, header=True)
		vcffiles[i] = annotatedVCFname

	#importantMuts = pd.DataFrame(listofmutstoExport)
	#print (importantMuts)
	#importantMuts = importantMuts.groupby(['segment', 'gene', 'position', 'AAchange'], as_index=False).mean()
	
	#importantMuts = importantMuts.loc[importantMuts['freq']>0.02 & importantMuts['freq']<0.98]
	#importantMuts['freq'] = importantMuts['freq']/len(vcffiles)
	#importantMutsexport = importantMuts.to_csv(os.path.join(outDir, "allMutationsPresent.tsv"), sep = '\t', index=None, header=True)


	#get end time
	end = time.time()
	#get total runtime
	runtimeSeconds = end - start1
	runtime = datetime.timedelta(seconds=runtimeSeconds)
	print(f'\nSniffles: Finished annotating snps in {str(runtime)}')
	return (vcffiles)