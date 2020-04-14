import subprocess as sub
import os
import shlex
import calldocker as cd
import multiprocessing as mp
import time, datetime
from sc import procTitle,checkexists

def removeDuplicates(runCFG,bam_files,threads='1'):
	#initial parameters
	outDir =runCFG['exec']['outdir']
	logfile = runCFG['exec']['logfile']
	checkexists(os.path.join(outDir,'rm_dups'))
	outDir = os.path.join(outDir,'rm_dups')

	#notify starting to remove duplicates
	procTitle('Remove Duplicates', runCFG)
	print('\nSniffles: Removing duplicate reads')
	#get time at start
	start = time.time()

	#generate commands
	cmds = []
	output_list = []
	for path in bam_files:
		full_path = os.path.abspath(path)
		file_name = os.path.basename(full_path)
		path = os.path.dirname(full_path)
		id = file_name.split(".")[0]

		#remove duplicate reads command
		cmd = f'java -Xmx2g -jar /tools/picard.jar MarkDuplicates I=/in_dir/{id}.bam O=/out_dir/{id}.bam REMOVE_DUPLICATES=true M=/out_dir/{id}.removeDupMetrics.txt'
		cmds.append(cmd)

		#add id to finished list
		output_list.append(os.path.join(outDir,f'{id}.bam'))

	#set up multiprocessing
	pool = mp.Pool(processes=threads)

	#denote start of remove duplicate reads in logs
	with open(logfile,'a') as outlog:
		outlog.write('***********\n')
		outlog.write('Removing Duplicates\n')
		#start multiprocessing
		results = pool.starmap_async(cd.call,[[cmd,'/reads',{path:"/in_dir",outDir:"/out_dir"}] for cmd in cmds])
		pool.close()
		pool.join()
		stdouts = results.get()
		for stdout in stdouts:
			outlog.write('-----------\n')
			outlog.write(stdout)
		#denote end of logs
		outlog.write('***********\n')
	#get time at end
	end = time.time()

	#determine runtime of processes
	runtime = round(end - start,2)
	print(f'\nSniffles: Finished removing duplicates in {runtime} seconds')
	return output_list


def normCoverage(runCFG,bam_files,threads='1'):
	#initial parameters
	outDir =runCFG['exec']['outdir']
	checkexists(os.path.join(outDir,'normalized'))
	logfile = runCFG['exec']['logfile']
	outDir = os.path.join(outDir,'normalized')

	#notify starting to remove duplicates
	procTitle("Downsampling with seqtk to normalize coverage", runCFG)
	#print('\n-----------------------Sniffles: Downsampling with seqtk to normalize coverage-----------------------')

	#get time at start
	start = time.time()

	#denote start of remove duplicate reads in logs
	with open(logfile,'a') as outlog:
		outlog.write('***********\n')
		outlog.write('Downsampling with seqtk to normalize coverage\n')

		#run normalization
		output_list = []
		for path in bam_files:
			
			full_path = os.path.abspath(path)
			file_name = os.path.basename(full_path)
			path = os.path.dirname(full_path)
			id = file_name.split('.')[0]

			#get reads from mapped bamfile
			cmd_get_reads = f'bash -c \'samtools collate /bam_files/{id}.bam collating && samtools fastq -n /bam_files/{id}.bam -1 /out_dir/{id}_mapped_1.fastq.gz -2 /out_dir/{id}_mapped_2.fastq.gz'

			#run seqtk to subsample readsc
			total_reads = runCFG['exec']['totalReads']
			cmd_normalization = f' && seqtk sample -s100 /out_dir/{id}_mapped_1.fastq.gz {total_reads} > {id}_1.fastq.gz && seqtk sample -s100 /out_dir/{id}_mapped_2.fastq.gz {total_reads} > /out_dir/{id}_2.fastq.gz\''
			
			if runCFG['exec']['unpaired']:
				cmd_get_reads += f' -0 /out_dir/{id}_mapped_U.fastq.gz && seqtk sample -s100 /out_dir/{id}_mapped_U.fastq.gz {total_reads} > /out_dir/{id}_U.fastq.gz'
				output_list.append([os.path.join(outDir,f'{id}_1.fastq.gz'),os.path.join(outDir,f'{id}_2.fastq.gz'),os.path.join(outDir,f'{id}_U.fastq.gz')])
			else:
				output_list.append([os.path.join(outDir,f'{id}_1.fastq.gz'),os.path.join(outDir,f'{id}_2.fastq.gz')])
			
			#start docker containers and run
			outlog.write(f'{id}\n-----------\n')
			stdout=cd.call(cmd_get_reads+cmd_normalization,'/out_dir',{path:"/bam_files",outDir:"/out_dir"})
			outlog.write(stdout)

			#cleanup
			try:
				os.remove(f'{outDir}/{id}_mapped_1.fastq.gz')
			except:
				pass
			try:
				os.remove(f'{outDir}/{id}_mapped_2.fastq.gz')
			except:
				pass
			try:
				os.remove(f'{outDir}/{id}_mapped_U.fastq.gz')
			except:
				pass
		outlog.write('***********\n')

	#get time at end
	end = time.time()

	#determine runtime of processes
	runtime = round(end - start,2)
	runtime = str(datetime.timedelta(seconds=runtime))
	print(f'\nSniffles: Finished normalizing read coverage in {runtime}')
	return output_list
