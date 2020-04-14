import os
import time, datetime
from sc import procTitle
import calldocker as cd
import multiprocessing as mp
import pandas as pd
from shutil import copytree, rmtree
import glob


def calltoVCF(callfile, outDir):
	#import .call file as pandas DF, extract sample name from .call name 
	callDF = pd.read_csv(callfile, sep="\t")
	samplename = "".join(callfile.split("/")[-1].split(".")[:-1])

	#create the mandatory VCF columns with universal values
	callDF["QUAL"] = "."
	callDF["INFO"] = "."
	callDF["FORMAT"] = "GT:DP:AD:FREQ"
	
	#create a "vcfDF" that consists of the .call SNPs that passed RePlow's quality metrics
	vcfDF = callDF.loc[callDF["Filter"] == "PASS"].reset_index()

	#Extract SNP depth, minor allele depth, allele frequency, 
	vcfDF["DP"] = vcfDF[["depth_mQ_rep1", "depth_mQ_rep2"]].apply(pd.to_numeric).sum(1).round(0).astype(int)
	vcfDF["FREQ"] = (vcfDF[["estimatedBAF_rep1", "estimatedBAF_rep2"]].apply(pd.to_numeric).mean(1))
	vcfDF["AD"] = (vcfDF["DP"]*vcfDF["FREQ"]).astype('int')#vcfDF[["bCnt_rep1", "bCnt_rep2"]].apply(pd.to_numeric).sum(1).round(0).astype(int)
	vcfDF["FREQ"] = (vcfDF["FREQ"]*100).round(4).astype(str) + "%"
	vcfDF["GT"] = "." #Genotype appears to be necessary for downstream VCFtools usage, but is not really applicable to our samples

	#reorder columns into appropriate final order
	vcfDF = vcfDF[["#chr", "pos", "ID", "ref", "alt", "QUAL", "Filter", "INFO", "FORMAT", "GT", "DP", "AD", "FREQ"]]

	#join GT, DP, AD, and FREQ columns into one sample column with those values separated by ":", as specified in 
	vcfDF[samplename] = vcfDF.loc[:,"GT":"FREQ"].astype(str).apply(":".join, axis=1)
	#Turn .call lowercase column names into VCF uppercase column names, and rearrange columns into final order
	vcfDF = vcfDF.rename(columns={"#chr":"#CHROM", "pos":"POS", "ref":"REF", "alt":"ALT", "Filter":"FILTER"})
	vcfDF = vcfDF[["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", samplename]]

	#export vcfDF as .vcf file with appropriate metadata to top of vcf
	#easiest way I found to do this was to export a TSV with just the metadata, then append vcfDF to bottom of that same file
	metadata = pd.Series(['##fileformat=VCFv4.2', '##source=JLL','##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\'Quality Read Depth of bases with Phred score >= 30\'>','##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\'Depth of variant-supporting bases (reads2)\'>', '##FORMAT=<ID=FREQ,Number=1,Type=String,Description=\'Variant allele frequency\'>'])

	outputvcfname = os.path.join(outDir, samplename + ".vcf")

	result = metadata.to_csv(outputvcfname, sep = "\t", index=False)

	banana = vcfDF.to_csv(outputvcfname, mode='a', sep="\t", index = False)
	
	return(samplename)

def RePlow(runCFG, bam_files, threads='1'):
	#set parameters
	outDir = runCFG['exec']['outdir']
	logfile = runCFG['exec']['logfile']
	bamfilespath = os.path.dirname(bam_files[0])#os.path.join(outDir, "norm_mapping")
	outDir = os.path.join(outDir,'snp_calls')

	reference_sequence_path = os.path.join(runCFG['exec']['outdir'], 'ref_sequence')
	reference_sequence_name = os.path.basename(runCFG['exec']['referenceSequence'])
	# if os.path.isdir(reference_sequence_path) and os.listdir("ref_sequence") != os.listdir(reference_sequence_path) :
	# 	rmtree (reference_sequence_path)
	# 	copytree("ref_sequence", reference_sequence_path)
	# if not os.path.isdir(reference_sequence_path):
	# 	copytree("ref_sequence", reference_sequence_path)
	#starting time point
	start =  time.time()
	procTitle('Analyzing SNPs with RePlow', runCFG)
	#print('\n-----------------------Sniffles: Calling SNPs with RePlow-----------------------')

	bams = []
	sample_list = []
	repDict = {}

	#create list of bam files to run
	print ("bam_files:")
	print (bam_files)
	
	for path in bam_files:
		full_path = os.path.abspath(path)
		file_name = os.path.basename(full_path)
		path = os.path.dirname(full_path)
		id = file_name.split(".")[0]
		sample_list.append(id)
		bams.append('/infile/'+file_name)
		print (file_name)

		#if processing replicate runs, create dictionary samplename:[list of replicate bam files for that sample]
		if runCFG['exec']['replicates']:
			repBreakdown = runCFG['exec']['replicateNotation'].split("_")
			repBreakdown = "_".join(repBreakdown[:-1])
			repBreakdown = repBreakdown[:repBreakdown.find(r"\d")]
			repBreakdown = repBreakdown.split("Sample")
			repKey = file_name[file_name.find(repBreakdown[0])+len(repBreakdown[0]):file_name.find(repBreakdown[1])]

			if repKey not in repDict.keys():
				repDict[repKey] = [file_name]
			else:
				repDict[repKey].append(file_name)
	print (repDict)


	#import SNP quality parameters from config
	snp_qual_threshold=runCFG['snpcalling']['snpQualityThreshold']
	consensus_frequency=runCFG['snpcalling']['consensusFrequency']
	mut_rate=runCFG['replow_settings']['mutrate']
	map_qual_threshold=runCFG['replow_settings']['mapquality']

	#make a default .bed file from reference that instructs replow to call SNPs on the whole genome
	chrom = []
	chromEnd = []
	with open(os.path.join(reference_sequence_path, reference_sequence_name),'r') as refseq:
		for line in refseq.readlines():
			if line[0] == ">":
				chrom.append(line[1:].rstrip())
				
			elif len(chrom) != len(chromEnd):
				chromEnd.append(len(line)-2)  #-1 to adjust for zero indexing in a bed file
			else:
				chromEnd[len(chrom)] =+ len(line)
		
		bedfile = pd.DataFrame()
		bedfile["chrom"] = chrom
		bedfile["chromStart"] = 0
		bedfile["chromEnd"] = pd.Series(chromEnd)
	print (reference_sequence_name)
	bedfilename = "".join(reference_sequence_name.split(".")[:-1])+".bed"
	bedcsv = bedfile.to_csv(os.path.join(reference_sequence_path, bedfilename), index=False, header=False, sep="\t")


	#Prep RePlow cmds
	cmds=[]
	keycmds=[]
	indexcmds=[]
	#Adds command to index reference sequence for RePlow
	indexcmds.append(f'samtools index /ref/{reference_sequence_name}')

	#add commands to list for multiprocessing
	###TO BE IMPROVED: Currently I only generate commands if I am running RePlow on replicate sequencing runs.
	### While that is the typical case, I should eventually expand this to cover single runs as well.
	#run through each sample name in replicate dictionary
	# print ("repDict:")
	print (repDict)
	for key in repDict.keys():
		#generate command to index each replicate bam file 
		for i, bam in enumerate(repDict[key]):
			repDict[key][i] = "/data/"+ bam
			indexcmds.append(f'samtools index {repDict[key][i]}')
		
		#create comma-deliniated list of replicate bam files, generate unique replow command for each sample
		bamslist=','.join(repDict[key])
		#replowcmd = f'java -cp dependency/*:classes tgil.replow.RePlow -r /ref/{reference_sequence_name} -b {bamslist} -T /ref/{bedfilename} -R /usr/bin/Rscript -f {consensus_frequency} -q {snp_qual_threshold} -Q {map_qual_threshold} -m {mut_rate} -o /output -L {key}'
		outlogHeader = f"\"{key}\n-----------\n\""
		replowcmd = f'bash -c \'printf {outlogHeader} >> {os.path.join("/logfile", os.path.basename(logfile))} && java -jar /source/RePlow-1.1.0.jar -r /ref/{reference_sequence_name} -b {bamslist} -T /ref/{bedfilename} -R /usr/bin/Rscript -f {consensus_frequency} -q {snp_qual_threshold} -Q {map_qual_threshold} -m {mut_rate} -o /output -L {key}\''
		keycmds.append(replowcmd)
	print ("Indexcmds:")
	print (indexcmds)
	print ("Keycmds:")
	print (keycmds)
	#generate mulitprocessing pool
	pool = mp.Pool(processes=threads)

	#index files first

	with open(logfile,'a') as outlog:
		outlog.write('***********\n')
		outlog.write('RePlow\n')
		outlog.write('***********\n')
		#run commands in docker contaniers with multiprocessing
		indexresults = pool.starmap_async(cd.call,[[cmd,'/source',{os.path.join (os.getcwd(), "replow"):"/source", bamfilespath:"/data", outDir:"/output", reference_sequence_path:"/ref",os.path.dirname(logfile):"/logfile"}] for cmd in indexcmds])
		pool.close()
		pool.join()

		pool = mp.Pool(processes=threads)
		results = pool.starmap_async(cd.call,[[cmd,'/source',{os.path.join (os.getcwd(), "replow"):"/source", bamfilespath:"/data", outDir:"/output", reference_sequence_path:"/ref",os.path.dirname(logfile):"/logfile"}] for cmd in keycmds])
		pool.close()
		pool.join()
		stdouts = indexresults.get() + results.get()
		
		for stdout in stdouts:
			#outlog.write('-----------\n')
			outlog.write(stdout)

		#Convert .call into .vcf
		for file in glob.glob(outDir + "/*.call"):
			samplename = calltoVCF(file, outDir)
			outlog.write(f"Call file {samplename} converted to VCF format")

		#denote end of logs
		outlog.write('***********\n')
				
	#get end time
	end = time.time()
	#get total runtime
	runtime = round(end - start,2)
	runtime = str(datetime.timedelta(seconds=runtime))
	print(f'\nSniffles: Finished calling snps in {runtime}')

	return (os.path.join(outDir, samplename +".vcf"))