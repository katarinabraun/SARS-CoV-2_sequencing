import os
import shlex
import subprocess as sub
import multiprocessing as mp
import time, datetime
from sc import procTitle
import calldocker as cd

def trim(readData,runCFG,threads=1,ids=''):

	#parameters
	minlength = runCFG['trimmomatic']['minlength']
	windowsize = runCFG['trimmomatic']['windowSize']
	qscore = runCFG['trimmomatic']['qscore']
	adapterpath = "/Trimmomatic-0.36/adapters/" + runCFG['trimmomatic']['adaptersFileName']
	outDir = runCFG['exec']['outdir']
	logfile = runCFG['exec']['logfile']

	#set up list of ids to trim
	if not ids:
		ids = readData.ids

	#generate commands for each trim job
	cmds = []
	#rmspacecmds = []
	for id in ids:
		#get read path
		if readData.reads[id]:
			read_path = os.path.dirname(os.path.abspath(readData.reads[id].fwd))
			read1_basename = os.path.basename(readData.reads[id].fwd)
			read2_basename = os.path.basename(readData.reads[id].rev)
		
		#main command
		outlogHeader = f"\"{id}\n-----------\n\""
		containerLogpath = os.path.join("/logfile", os.path.basename(logfile))
		main_cmd = f'bash -c \'printf {outlogHeader} >> {containerLogpath} && java -jar /tools/trimmomatic.jar '
		#regexExpression='s/(^@.*) (.*)/\\1_\\2/g'

		#determine args
		if runCFG['trimmomatic']['removeAdapters']:
			if runCFG['trimmomatic']['paired']:
				args = f'PE {read1_basename} {read2_basename} -baseout /output/{id}_trimmed.fastq.gz ILLUMINACLIP:{adapterpath}:1:30:10 SLIDINGWINDOW:{windowsize}:{qscore} MINLEN:{minlength}\''
				readData.add_runtime('trimmed',id,f'{outDir}/trimmed/{id}_trimmed_1P.fastq.gz',f'{outDir}/trimmed/{id}_trimmed_2P.fastq.gz')
				#rmspaces = f"bash -c \'sed -re \"{regexExpression}\" -i /output/{id}_trimmed_1P.fastq.gz && sed -re \"{regexExpression}\" -i /output/{id}_trimmed_2P.fastq\'"
			else:
				args = f'SE {read1_basename} -baseout /output/{id}_trimmed.fastq.gz ILLUMINACLIP:{adapterpath}:1:30:10 SLIDINGWINDOW:{windowsize}:{qscore} MINLEN:{minlength}\''
				readData.add_runtime('trimmed',id,f'{outDir}/trimmed/{id}_trimmed.fastq.gz')
				#rmspaces = f"bash -c \'sed -re \"{regexExpression}\" -i /output/{id}_trimmed.fastq\'"
		else:
			if runCFG['trimmomatic']['paired']:
				args = f'PE {read1_basename} {read2_basename} -baseout /output/{id}_trimmed.fastq.gz SLIDINGWINDOW:{windowsize}:{qscore} MINLEN:{minlength}\''
				readData.add_runtime('trimmed',id,f'{outDir}/trimmed/{id}_trimmed_1P.fastq.gz',f'{outDir}/trimmed/{id}_trimmed_2P.fastq.gz')
				#rmspaces = f"bash -c \'sed -re \"{regexExpression}\" -i /output/{id}_trimmed_1P.fastq && sed -re \"{regexExpression}\" -i /output/{id}_trimmed_2P.fastq\'"
			else:
				args = f'SE {read1_basename} -baseout /output/{id}_trimmed.fastq.gz SLIDINGWINDOW:{windowsize}:{qscore} MINLEN:{minlength}\''
				readData.add_runtime('trimmed',id,f'{outDir}/trimmed/{id}_trimmed.fastq.gz')
				#rmspaces = f"bash -c \'sed -re \"{regexExpression}\" -i /output/{id}_trimmed.fastq\'"

		#prepare command and add to list
		sample_cmd = main_cmd+args
		#rmspacecmds.append(rmspaces)
		cmds.append(sample_cmd)
		for cmd in cmds:
			print (cmd)

	#make out dir if it doesn't already exist
	try:
		os.mkdir(os.path.join(outDir,'trimmed'))
	except:
		pass

	#set up multiprocessing
	pool = mp.Pool(processes=threads)
	# pool2 = mp.Pool(processes=threads)
	#notify starting trimming
	procTitle("Started quality trimming", runCFG)
	
	#start timer
	start = time.time()

	#denote logs
	with open(logfile,'a') as outlog:
		outlog.write('***********\n')
		outlog.write('Trimmomatic\n')
		outlog.write('***********\n')
		#begin multiprocessing
		results = pool.starmap_async(cd.call,[[cmd,'/data',{read_path:"/data",os.path.join(outDir,'trimmed'):"/output",os.path.dirname(logfile):"/logfile"}] for cmd in cmds])
		pool.close()
		pool.join()
		stdouts = results.get()
		for stdout in stdouts:
			outlog.write(stdout)
		# results = pool2.starmap_async(cd.call,[[cmd,'/data',{read_path:"/data",os.path.join(outDir,'trimmed'):"/output",os.path.dirname(logfile):"/logfile"}] for cmd in rmspacecmds])
		# pool.close()
		# pool.join()
		# stdouts = results.get()
		# for stdout in stdouts:
		# 	outlog.write(stdout)
		#denote end of logs
		outlog.write('***********\n')

	#get time
	end = time.time()
	#determine runtime of processes
	runtime = round(end - start,2)
	runtime = str(datetime.timedelta(seconds=runtime))
	print(f'Sniffles: Finished trimming in {runtime} seconds\n')
