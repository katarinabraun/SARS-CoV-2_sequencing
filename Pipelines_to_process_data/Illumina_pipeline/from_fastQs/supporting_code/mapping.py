import os, re
import multiprocessing as mp
import time, datetime
from sc import procTitle,checkexists,cpu_count
import calldocker as cd
from shutil import copyfile


seed = 0

def indexing(runCFG, *paths):
	#print('\n-----------------------Sniffles: Indexing reference sequence-----------------------\n')

	logfile = runCFG['exec']['logfile']
	outDir = runCFG['exec']['outdir'] + '/ref_sequence'
	checkexists(outDir)
	procTitle("Indexing reference sequence", runCFG)

	for path in paths:
		reference_sequence_abspath = os.path.abspath(path)
		reference_sequence_name = os.path.basename(reference_sequence_abspath)
		print (path)
		#index reference
		cmd = f'bowtie2-build {reference_sequence_name} {reference_sequence_name} --quiet'#' --threads {threads}'
		with open(logfile,'a') as outlog:
			outlog.write("***********\n")
			outlog.write("Bowtie2 indexing the reference\n")
			copyfile(reference_sequence_abspath,os.path.join(outDir,reference_sequence_name))
			outlog.write(cd.call(cmd,'/data',{outDir:"/data"}))
			outlog.write("***********\n")

#mapping parameters
#param_path - list of paths in the following order
#   [(id,path to fwd read, path to rev read, path to reference),...]
#   if not paired end then "path to rev read" will be empty
#threads - number of threads
#outDir - the output directory

def mapping(runCFG,param_paths,outDir,threads='8'):
	
	seed = 0
	start=time.time()

	#assert len(param_paths) > 0, "Cannot map reads: No reads provided"
	print (param_paths)
	logfile = runCFG['exec']['logfile']

	num_jobs,num_threads = cpu_count(threads)
	cmds = []
	read_path = ''
	ref_path = ''
	output_bam_list = []
	for param_path in param_paths:
		id = param_path[0]
		read1 = os.path.basename(param_path[1])
		read2 = os.path.basename(param_path[2])
		try:
			unpairedReads = os.path.basename(param_path[4])
		except:
			unpairedReads = ""
		print (f"unpairedReads for id {id}: " + unpairedReads)
		read1un = read1.split(".")[0][0:-1]+"U.fastq.gz"
		read2un = read2.split(".")[0][0:-1]+"U.fastq.gz"
		read_path = os.path.dirname(os.path.abspath(param_path[1]))
		ref_path = runCFG['exec']['outdir'] + '/ref_sequence'
		reference_sequence_name = os.path.basename(param_path[3])
		outlogHeader = f"\"{id}\n-----------\n\""
		#check output folder exists
		checkexists(os.path.join(outDir))
		logfilepath = os.path.join("/logfile", os.path.basename(logfile))
		checkexists(os.path.join(outDir, "unmapped"))

		
		if read2 != '':
			#generate command for paired end
			cmd = f"bash -c \'printf {outlogHeader} >> {logfilepath} && bowtie2 -x {reference_sequence_name} -1 /reads/{read1} -2 /reads/{read2} --seed {seed} --no-unal --local 2>> /logfile/{id}mappingstats.log | samtools view -bS | samtools sort -o /output/{id}"
		else:
			#generate command for interleaved
			cmd = f"bash -c \'printf {outlogHeader} >> {logfilepath} && bowtie2 -x {reference_sequence_name} --interleaved /reads/{read1} --seed {seed} --no-unal --local 2>> /logfile/{id}mappingstats.log | samtools view -bS | samtools sort -o /output/{id}"
		
		if runCFG['exec']['unpaired']:
			if unpairedReads != "":
				cmdUn = f"_paired.bam && printf {outlogHeader} >> {logfilepath} && bowtie2 -x {reference_sequence_name} -U /reads/{unpairedReads} --un /output/unmapped/U_{read1} --seed {seed} --no-unal --local 2>> /logfile/{id}mappingstats.log | samtools view -bS | samtools sort -o /output/{id}_unpaired.bam"
			elif read2 != '':
				#generate command for unpaired
				cmdUn = f"_paired.bam && printf {outlogHeader} >> {logfilepath} && bowtie2 -x {reference_sequence_name} -U /reads/{read1un},/reads/{read2un} --un /output/unmapped/U_{read1} --seed {seed} --no-unal --local 2>> /logfile/{id}mappingstats.log | samtools view -bS | samtools sort -o /output/{id}_unpaired.bam"
			else:
				#generate command for interleaved
				cmdUn = f"_paired.bam && printf {outlogHeader} >> {logfilepath} && bowtie2 -x {reference_sequence_name} -U /reads/{read1un} --un /output/unmapped/U_{read1} --seed {seed} --no-unal --local 2>> /logfile/{id}mappingstats.log | samtools view -bS | samtools sort -o /output/{id}_unpaired.bam"
		
			mergecmd = f" && samtools merge /output/{id}.bam /output/{id}_unpaired.bam /output/{id}_paired.bam && samtools index /output/{id}.bam\'"
			cmds.append(cmd+cmdUn+mergecmd)
		else:
			cmds.append(cmd+".bam\'")
		#data for next stage
		output_bam_list.append(os.path.join(outDir,f'{id}.bam'))

	#set up multiprocessing
	pool = mp.Pool(processes=num_jobs)

	#get start time
	start = time.time()
	
	#denote start of mapping in logs
	with open(logfile,'a') as outlog:
		outlog.write('***********\n')
		outlog.write('Mapping\n')
		#start multiprocessing
		results = pool.starmap_async(cd.call,[[cmd,'/reads',{ref_path:"/reference",read_path:"/reads",outDir:"/output",os.path.dirname(logfile):"/logfile"}] for cmd in cmds])
		pool.close()
		pool.join()
		stdouts = results.get()
		for stdout in stdouts:
			outlog.write('-----------\n')
			outlog.write(stdout)
		#denote end of logs
		outlog.write('***********\n')
	
	#get end time
	end = time.time()
	#get total runtime
	runtime = round(end - start,2)
	runtime = str(datetime.timedelta(seconds=runtime))
	print(f'\nSniffles: Finished mapping in {runtime}')
	return output_bam_list

def average_depth(runCFG,bam_list,inDir,outDir):
	print('\n************* Filtering reads by minimum depth ************')

	#check that output folder exists
	checkexists(os.path.join(outDir))

	#setup initial parameters
	ref_path = runCFG['exec']['outdir'] + '/ref_sequence'
	reference_sequence_name = os.path.basename(runCFG['exec']['referenceSequence'])
	logfile = runCFG['exec']['logfile']

	#bam file list that will meet threshold
	filtered_bam_list = []

	#create logfile that will hold all average depths
	with open(os.path.join(outDir,'average_depth.log'),'w') as outdepth:

		#loop through each bam
		for bam in bam_list:
			filename = os.path.basename(bam)
			id = filename.split('.')[0]
			#open the log to log output
			with open(logfile,'a') as outlog:
				outlog.write('***********\n')
				outlog.write('Determining Average Coverage\n')
				outlog.write('***********\n')
				outlog.write(f'{id}\n-----------\n')
				#generate command for the current bam file
				cmd = f'bash -c \'samtools view /indata/{filename} > {id}.tmp.sam  && /tools/bbmap/pileup.sh in={id}.tmp.sam out={id}_coverage.tsv ref=/reference/{reference_sequence_name} && rm {id}.tmp.sam\''
				#use docker to run the command
				output = cd.call(cmd,'/outdata',{inDir:"/indata",outDir:"/outdata",ref_path:"/reference"})
				#record the output in the log
				outlog.write(output)
				#denote end of logs
				outlog.write('***********\n')

			#only add isolates that pass average depth and percent of reference covered
			percent_cov = 0
			avg_cov = 0
			#parse lines of stdout for info we need
			for line in output.splitlines():
				if "Percent of reference bases covered:" in line:
					match = re.search('[0-9,.]+',line.split(":")[-1])
					if match:
						percent_cov = float(match[0])

				if "Average coverage:" in line:
					match = re.search('[0-9,.]+',line.split(":")[-1])
					if match:
						avg_cov = float(match[0])
			#check against the config for min thresholds
			if avg_cov >= runCFG['exec']['minimumAverageDepth'] and percent_cov >= runCFG['exec']['percentRefCovered']:
				#record
				outdepth.write(f'{id},{avg_cov},{percent_cov},Pass\n')
				filtered_bam_list.append(bam)
			else:
				outdepth.write(f'{id},{avg_cov},{percent_cov},Fail\n')
				print (f'{id},{avg_cov},{percent_cov},Fail\n')
					
	#return only bam files that meet coverage requirements
	return filtered_bam_list
	