import os
import shlex
import subprocess as sub
import time
from sc import procTitle,checkexists
import multiprocessing as mp
import calldocker as cd
from Bio import SeqIO
from Bio.Seq import Seq
import random
from collections import defaultdict
from operator import itemgetter
from subprocess import run
import glob

ambNuc = {"N":("A","T","G","C"), 'W':('A','T'), 'S':('G', 'C'), 'B':('T','G','C'), 'D':('A','T','G'), 'H':('A','T','C'), 'V':('A','G','C'), 'K':('G','T'), 'M':('A','C'), 'R':('A','G'), 'Y':('C','T')}

def consensus(runCFG,bam_list,threads='1'):
	#fastas have been initially mapped, then checked for coverage
	#then downsampled, then mapped again
	#now we will generate a consensus to map against

	#initial parameters
	outDir =runCFG['exec']['outdir']
	logfile = runCFG['exec']['logfile']
	outDir = os.path.join(outDir,'consensus')
	checkexists(outDir)

	#get start time
	overall_start = time.time()
	start = time.time()

	procTitle("Generating consensus sequence", runCFG)
	

	#set reference sequence
	reference_sequence_abspath = os.path.abspath(runCFG['exec']['referenceSequence'])
	reference_sequence_name = os.path.basename(reference_sequence_abspath)
	reference_sequence_dir = runCFG['exec']['outdir'] + '/ref_sequence'
	reference_sequence_id = reference_sequence_name.split(".")[0]

	#command list
	cmds = []
	vcf_list = []
	for path in bam_list:
		full_path = os.path.abspath(path)
		file_name = os.path.basename(full_path)
		path = os.path.dirname(full_path)
		id = file_name.split(".")[0]
		minCov = runCFG['snpcalling']['minCoverage']
		quality = runCFG['snpcalling']['snpQualityThreshold']
		freq = runCFG['snpcalling']['consensusFrequency']

		#make multiway pileup using samtools
		cmd1 = f'bash -c \'samtools mpileup -ABd 1000000 /infile/{file_name} -f /ref/{reference_sequence_name} -o /outfile/{id}.pileup'
		
		#run varscan mpileup2cns to generate vcf with consensus information
		cmd2 = f'java -jar /tools/varscan.jar mpileup2cns {id}.pileup --min-coverage {minCov} --min-avg-qual {quality} --min-var-freq {freq} --strand-filter 1 --variants --output-vcf 1 > {id}.vcf'
		
		cmds.append(cmd1 + " && " + cmd2 +"\'")
		vcf_list.append(os.path.join(outDir,f'{id}.vcf'))

	#setup multiprocessing
	pool = mp.Pool(processes=threads)
	#start multiprocessing
	with open(logfile,'a') as outlog:
		outlog.write('***********\n')
		outlog.write('1st round of SNP calling\n')
		#start multiprocessing
		results = pool.starmap_async(cd.call,[[cmd,'/outfile',{reference_sequence_dir:"/ref",path:"/infile",outDir:"/outfile"}] for cmd in cmds])
		pool.close()
		pool.join()
		stdouts = results.get()
		for stdout in stdouts:
			outlog.write('-----------\n')
			outlog.write(stdout)
		#denote end of logs
		outlog.write('***********\n')

	#check if vcf file is empty, if it is skip id and remove vcf file
	filtered_vcf_list = []
	for path in vcf_list:
		try:
			if os.path.getsize(path)>0:
				filtered_vcf_list.append(path)
			else:
				os.remove(path)
		except:
			pass

	end = time.time()
	runtime = round(end - start,2)
	print(f'\nSniffles: Finished generating the consensus vcf in {runtime} seconds\n')
	start = time.time()

	print(f'\nSniffles: Generating consensus fasta\n')
	#command list for compressing files
	cmds = []
	out_fasta = []
	merge_fastas = []
	for vcf in filtered_vcf_list:
		full_path = os.path.abspath(vcf)
		file_name = os.path.basename(full_path)
		path = os.path.dirname(full_path)
		id = file_name.split(".")[0]
		if len(bam_list) == 1:
			id = id+"_consensus_noambig"
		#compress vcf file with bgzip and create consensus for each sample
		cmd = f'bash -c \'bgzip {id}.vcf && tabix {id}.vcf.gz && bcftools norm -Ob -m -any -f /ref/{reference_sequence_name} {id}.vcf.gz -o {id}_norm.bcf.gz && tabix {id}_norm.bcf.gz && bcftools consensus -f /ref/{reference_sequence_name} {id}_norm.bcf.gz -o {id}.fasta\''
		out_fasta.append(os.path.join(outDir,f'{id}.fasta'))
		merge_fastas.append(f'{id}.fasta')
		cmds.append(cmd)

	pool = mp.Pool(processes=threads)
	#create consensus for all samples in reference
	with open(logfile,'a') as outlog:
		outlog.write('***********\n')
		outlog.write('Creating individual consensus fastas\n')
		#start multiprocessing
		results = pool.starmap_async(cd.call,[[cmd,'/outfile',{reference_sequence_dir:"/ref",outDir:"/outfile"}] for cmd in cmds])
		stdouts = results.get()
		pool.close()
		pool.join()
		for stdout in stdouts:
			outlog.write('-----------\n')
			outlog.write(stdout)
		#denote end of logs
		outlog.write('***********\n')
	
	if len(bam_list) == 1:
		consensusRef = id+'.fasta'
		return consensusRef

	#then separate by segment and create consensus for each segment
	consensuscmds = []
	segfiles=[]
	consensusRef = (reference_sequence_id)+"_consensus.fasta"
	sequencelist=defaultdict(list)
	for fastaseq in merge_fastas:
		with open(os.path.join(outDir, fastaseq), "r") as fs:
			for segment in SeqIO.parse(fs, "fasta"):
				segname = segment.id.replace('/','_')
				segment.id = fastaseq.split("/")[-1].split(".")[0].split("_")[0]
				segment.name = segname.split("_")[-1]
				segment.description = ""
				sequencelist[segname].append(segment)
	for seg in sequencelist.keys():
		segmentfile = seg+".fasta"
		segfile = segmentfile.split(".")[0]+"_consensus.fasta"
		segfiles.append(segfile)
		with open(os.path.join(outDir, segmentfile), "w") as f:
			SeqIO.write(sequencelist[seg], f, "fasta-2line")
		consensuscmds.append(f'bash -c \'clustalo --infile {segmentfile} | consambig -filter -name {segmentfile.split(".")[0]} > {segfile}\'')
		
	pool = mp.Pool(processes=threads)
	#start multiprocessing
	with open(logfile,'a') as outlog:
		outlog.write('***********\n')
		outlog.write('Creating population consensus fasta\n')
		#start multiprocessing
		results = pool.starmap_async(cd.call,[[cmd,'/outfile',{reference_sequence_dir:"/ref",outDir:"/outfile"}] for cmd in consensuscmds])
		pool.close()
		pool.join()
		stdouts = results.get()
		for stdout in stdouts:
			outlog.write('-----------\n')
			outlog.write(stdout)
		#denote end of logs
		outlog.write('***********\n')

	#concatenate all segment consensuses into one reference consensus
	with open (os.path.join(outDir, consensusRef), 'a+') as confile:
		for segfile in segfiles:
			with open(os.path.join(outDir, segfile), 'r') as segmentfasta:
				confile.write(segmentfasta.read())


	#correct ambiguous bases
	outputConRefs = []
	with open(os.path.join(outDir, consensusRef), "r") as ambig:
		refreads = SeqIO.parse(ambig, 'fasta')
		for segment in refreads:
			segSeq = list(segment)
			for nuc in ambNuc.keys():
				if nuc in segSeq:
					print (f"Degenerate nucleotide {nuc} is in {segment.id} at locations {[i for i, x in enumerate(segSeq) if x == nuc]}")
					degNucs = [i for i, x in enumerate(segSeq) if x == nuc]
					with open (os.path.join(outDir, segment.id+".fasta")) as sampleseq:
						sampleSeqs = list(SeqIO.parse(sampleseq, 'fasta'))
						for i, ix in enumerate(degNucs):
							degNucs[i] = (ix, [seq[ix] for seq in sampleSeqs])
						for nuc in degNucs:
							#determine most common nucleotide; if more than one most common, choose at random
							b = defaultdict(int)
							d = defaultdict(list)
							for key, value in enumerate(nuc[1]):
								b[value] += 1
							for key, value in b.items():
								d[value].append(key)
							segSeq[nuc[0]] = random.choice(max(d.items())[1]).upper()
			segment.seq=Seq("".join(segSeq))
			outputConRefs.append(segment)
	consensusRef = os.path.join(outDir, consensusRef.split(".")[0]+"_noambig.fasta")
	runCFG['exec']['referenceSequence'] = consensusRef
	with open(consensusRef, "w") as conRef:
		SeqIO.write(outputConRefs, conRef, 'fasta-2line')

	end = time.time()
	runtime = round(end - start,2)
	print(f'\nSniffles: Finished generating consensus fasta in {runtime} seconds')

	#determine runtime of processes
	end = time.time()
	runtime = round(end - overall_start,2)
	print(f'\nSniffles: Finished generating consensus sequence in {runtime} seconds')
	return consensusRef