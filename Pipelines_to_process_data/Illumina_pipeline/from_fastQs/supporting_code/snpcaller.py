import os
import shlex
import subprocess as sub
import multiprocessing as mp
import time, datetime
from sc import procTitle,checkexists
import calldocker as cd
from shutil import copytree, rmtree
from VCFaverager import VCFaverager
from VCFcombiner import VCFcombiner

def snpcaller(runCFG,bam_files,threads='1'):
	#set parameters
	outDir = runCFG['exec']['outdir']
	logfile = runCFG['exec']['logfile']
	outDir = os.path.join(outDir,'snp_calls')
	checkexists(outDir)
	
	#set reference sequence
	reference_sequence_path = os.path.dirname(runCFG['exec']['referenceSequence'])
	reference_sequence_name = os.path.basename(runCFG['exec']['referenceSequence'])
	

	#starting time point
	start =  time.time()
	if runCFG['exec']['replicates']:
		message = 'Calling replicate SNPs with Varscan'
	else:
		message = 'Calling SNPs with Varscan'

	procTitle(message, runCFG)

	bams = []
	sample_list = []
	listofVCFs = []
	repDict = {}
	
	#Create list of bam files to call SNPs on
	for path in bam_files:
		full_path = os.path.abspath(path)
		file_name = os.path.basename(full_path)
		path = os.path.dirname(full_path)
		id = file_name.split(".")[0]
		sample_list.append(id)
		bams.append('/infile/'+file_name)

		#if processing replicate runs, create dictionary samplename:[list of replicate vcf files for that sample]
		#this dictionary will be used later to merge and average replicate vcfs
		if runCFG['exec']['replicates']:
			repBreakdown = runCFG['exec']['replicateNotation'].split("_")
			repBreakdown = "_".join(repBreakdown[:-1])
			repBreakdown = repBreakdown[:repBreakdown.find(r"\d")]
			repBreakdown = repBreakdown.split("Sample")
			repKey = file_name[file_name.find(repBreakdown[0])+len(repBreakdown[0]):file_name.find(repBreakdown[1])]

			vcf_name = (id+".vcf")
			listofVCFs.append(vcf_name)
			if repKey not in repDict.keys():
				repDict[repKey] = [vcf_name]
			else:
				repDict[repKey].append(vcf_name)

	#import SNP calling quality parameters from config file
	snp_frequency=runCFG['snpcalling']['snpFrequency']
	min_cov=runCFG['snpcalling']['minCoverage']
	snp_qual_threshold=runCFG['snpcalling']['snpQualityThreshold']

	#generate commands to call variants
	cmds=[]
	for bam, sample in zip(bams, sample_list):
		#mpileup command
		outlogHeader = f"{bam.split('/')[-1].split('.')[0]}\n-----------\n"

		cmd1 = f'printf \"{outlogHeader}\" >> {os.path.join("/logfile", os.path.basename(logfile))} && samtools mpileup -ABR -d 1000000 {bam} -f /ref/{reference_sequence_name} > {sample}.mpileup'

		#varscan command
		cmd2 = f'java -jar /tools/varscan.jar mpileup2snp {sample}.mpileup --min-coverage {min_cov} --min-avg-qual {snp_qual_threshold} --min-var-freq {snp_frequency} --strand-filter 1 --output-vcf 1 --variants --vcf-sample-list <(echo -e "{sample}") > {sample}_temp.vcf'
		#compress and normalize vcf
		cmd3 = f'bcftools norm -c sw -m - -f /ref/{reference_sequence_name} -o {sample}.vcf {sample}_temp.vcf && rm {sample}_temp.vcf'
		
		if not runCFG['exec']['replicates']:
			listofVCFs.append(os.path.join(outDir, f"{sample}.vcf"))
		#add commands to list for multiprocessing
		cmds.append("bash -c \'" + cmd1 + " && " + cmd2 + " && " + cmd3 + "\'")
	
	#initialize multiprocessing pool
	pool = mp.Pool(processes=threads)

	#open logfile
	with open(logfile,'a') as outlog:
		outlog.write('***********\n')
		outlog.write('Calling SNPs\n')

		#run commands with mutliprocessing
		results = pool.starmap_async(cd.call,[[cmd, '/outfile',{reference_sequence_path:"/ref",path:"/infile",outDir:"/outfile",os.path.dirname(logfile):"/logfile"}] for cmd in cmds])
		
		pool.close()
		pool.join()
		stdouts = results.get()
		print ('finished all results')
		for stdout in stdouts:
			#outlog.write('-----------\n')
			outlog.write(stdout)
		#if processing duplicate runs, merge and average SNP calls
		if runCFG['exec']['replicates']:
			listofVCFs = (VCFaverager(runCFG, repDict, listofVCFs))

		outlog.write(str(listofVCFs))
		outlog.write('-----------\n')
		
		#Combine sample vcfs into one master VCF:
		#allSNPs = VCFcombiner(runCFG, listofVCFs, "allVarscanSNVs.vcf")
		#outlog.write(f"\nCombined all vcf files into master vcf file allVarscanSNVs.vcf\n")
		outlog.write('-----------\n')
		#denote end of logs
		outlog.write('***********\n')

	#get end time
	end = time.time()
	#get total runtime
	runtime = round(end - start,2)
	runtime = str(datetime.timedelta(seconds=runtime))
	print(f'\nSniffles: Finished calling snps in {runtime}')
	
	return (listofVCFs)
