#author: Kelsey Florek
#email: kelsey.florek@slh.wisc.edu
#stripped down docker calling function

import docker
import os
import time

daystring = "%a, %d %b %Y"
timestring = "%H:%M:%S"

def call(command,cwd='',paths={},container='jlalli/sniffles_plus',remove=True):
	###access docker environment
	client = docker.from_env()
	print (command)
	###get the effective user and group id's
	user = str(os.geteuid())+':'+str(os.getegid())
	
	###setup mount point paths
	#{"/path/outside":"/path/incontainer"}
	volumes = {}
	if paths:
		for key in paths.keys():
			volumes[key] = {'bind':paths[key],'mode':'rw'}
	print (volumes)
	###run the container
	#create empty variable for holding byte object output for the container logs
	output = b''
	
	#try block to run the container
	try:
		container_obj = client.containers.run(container, command, user=user,volumes=volumes,working_dir=cwd,remove=remove,detach=True,stderr=True,environment=["BOWTIE2_INDEXES=/reference/"])
	except:
		#loop through output as it is streamed
		for line in container_obj.logs(stream=True):
			output += time.strftime(timestring, time.localtime()).encode('ascii') + "\t".encode('ascii') + line
			print (line.decode('utf-8'))
			
	else:
		for line in container_obj.logs(stream=True):
			output += time.strftime(timestring, time.localtime()).encode('ascii') + "\t".encode('ascii') + line
			print (line.decode('utf-8'))
	#once container is finished return output as a string
	return output.decode('utf-8')
