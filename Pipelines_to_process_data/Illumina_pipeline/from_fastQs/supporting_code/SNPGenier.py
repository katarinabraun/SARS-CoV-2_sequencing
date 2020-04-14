from __future__ import division
import os
import subprocess as sub
from subprocess import run, PIPE
import multiprocessing as mp
import time, datetime
#from sc import procTitle
import calldocker as cd
import numpy as np
import shutil
import glob
from Bio import SeqIO
from VCFcombiner import VCFcombiner



def SNPgenier(runCFG, vcffolder, slidingwindow=None, threads=mp.cpu_count()):
	#set parameters
	outDir = os.path.join(runCFG['exec']['outdir'], 'SNPGenie_output')
	logfile = runCFG['exec']['logfile']

	if slidingwindow == None:
		  #sliding window size for SNPGenie analysis (# of codons)
		  slidingwindow = runCFG['postprocessing']['slidingwindow']

	#set reference sequence
	if "/" not in runCFG['exec']['referenceSequence']:
		reference_sequence_path = os.path.join(runCFG['exec']['outdir'],"ref_sequence", runCFG['exec']['referenceSequence'])
	else:
		reference_sequence_path = runCFG['exec']['referenceSequence']
	

	#starting time point
	start2 =  time.time()
	#procTitle("Performing population analysis with SNPGenie", runCFG)
	
	#delete old results folder from any prior run that was aborted
	if os.path.isdir("SNPGenie_Results"):
		print ("Found old SNPGenie results folder")
		shutil.rmtree("SNPGenie_Results")

	#generate individual segment fasta files and run snpgenie
	cmds = []
	segOutDirs=[]
	
	
	#use vcf-tools to merge all VCFs into one
	vcffiles = glob.glob(vcffolder + "/*1.vcf")
	
	if len(glob.glob(vcffolder+"/all*.vcf"))!=1:
		master_vcf_name = VCFcombiner(runCFG, vcffiles, "all_snps.vcf")
	else:
		master_vcf_name = glob.glob(vcffolder+"/all*.vcf")[0]
		
	cmds = []
	#separate out gtf, vcf, and reference fasta into individual segments in their own folders
	with open(reference_sequence_path, 'r') as refseqfile:
		refseqs = SeqIO.parse(refseqfile, "fasta")
		
		with open(runCFG['postprocessing']['gtfFileName'], 'r') as master_gtf:
			master_gtf_lines = master_gtf.read().splitlines()
			
			with open (os.path.join(vcffolder, master_vcf_name), 'r') as master_vcf:
				master_vcf_lines =  master_vcf.read().splitlines()

				for fastaseq in refseqs:
					segname = fastaseq.id.replace("/", "_")
					with open(segname+".fasta", "w") as f:
						SeqIO.write(fastaseq, f, "fasta-2line")
					with open(segname+".gtf", "w") as gtf:
						for line in master_gtf_lines:
							line = line.replace("/", "_")
							if segname in line:
								gtf.write(line+"\n")
					with open(segname+".vcf", "w") as vcfNew:
						vcfNew.write(master_vcf_lines[0]+"\n")
						for line in master_vcf_lines:
							line = line.replace("/", "_")
							if "#" in line:
								vcfNew.write(line+"\n")
							elif segname in line:
								vcfNew.write(line+"\n")
					
					segmentdir = "/segmentdir"
					
					segOutDir = os.path.join(outDir, segname)
					print (segOutDir)
					if not os.path.isdir(segOutDir):
						os.makedirs(segOutDir)
					
					shutil.move(segname+".fasta",os.path.join(outDir, segname, segname + ".fasta"))
					shutil.move(segname+".gtf",os.path.join(outDir, segname, segname + ".gtf"))
					shutil.move(segname+".vcf",os.path.join(outDir, segname, segname + ".vcf"))

					#copy snpgenie.pl into segment folder, allowing me to run snpgenie from that folder which makes setting directories much easier.
					shutil.copyfile("snpgenie.pl", os.path.join(segOutDir, "snpgenie.pl"))

					#create commands to run snpgenie on each segment
					outlogHeader = f"\"{segname}\n-----------\n\""
					cmd = (" ".join(("perl", "snpgenie.pl", "--vcfformat=4", "--snpreport=" + segmentdir +"/" + segname + ".vcf", "--fastafile="+ segmentdir +"/" +  segname + ".fasta", "--gtffile="+ segmentdir +"/" + segname + ".gtf", "--slidingwindow="+str(slidingwindow), "--minfreq="+str(runCFG['postprocessing']['minSNPfreq']), '--outputfolder='+segOutDir+"\'")))
					cmd = f'bash -c \'printf {outlogHeader} >> {os.path.join("/logfile", os.path.basename(logfile))} && ' +cmd
					cmds.append(cmd)
					segOutDirs.append(segOutDir)

	
	cmds = tuple(cmds)
	segOutDirs = tuple(segOutDirs)
	print (segOutDirs)
	#set up mutltiprocessing pool
	pool = mp.Pool(processes=threads)
	
	callcmds = []
	#create docker commands that have segment-specific directories
	for command, segpath in zip(cmds, segOutDirs):
		callcmds.append((command,'/segmentdir', {segpath:"/segmentdir",os.path.dirname(logfile):"/logfile"}, "acdaic4v/ubuntu-perl-base"))
	
	with open(logfile,'a') as outlog:
		outlog.write('***********\n')
		outlog.write('Getting Population Statistics\n')
		results = pool.starmap_async(cd.call, [callcmds[i] for i in range(0, len(callcmds))])
		pool.close()
		pool.join()
		try:
			stdouts = results.get()
			for stdout in stdouts:
				outlog.write('-----------\n')
				outlog.write(stdout)
			#denote end of logs
		except:
			pass
		outlog.write('***********\n')

	#get end time
	end = time.time()
	#get total runtime
	runtime = round(end - start2,2)
	runtime = str(datetime.timedelta(seconds=runtime))
	print(f'\nSniffles: Finished Getting Population Statistics in {runtime}')
	return (segOutDirs)