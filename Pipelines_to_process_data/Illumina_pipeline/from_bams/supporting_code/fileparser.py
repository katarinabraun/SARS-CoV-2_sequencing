import os

#class for containing all sequencing run information
class RunFiles:
	#class containing all of the read information for an isolate
	class Fastqs:
		#isolate id
		id = ''
		#path to forward reads fastq
		fwd = ''
		#path to reverse reads fastq
		rev = ''
		#path to reads fastq if unpaired or interleaved
		path = ''
		#switch for interleaved
		inteleaved = False
		#switch for paired or unpaired
		paired = True
		def __init__(self,id='',fwd='',rev='',path='',interleaved=False):
			#set all variables to definition arguments
			self.id = id
			self.fwd = fwd
			self.rev = rev
			self.path = path
			self.interleaved = interleaved
			if path:
				self.paired = False

	#data dictonary containing all of the read objects
	reads = {}

	#data dictonary for containing runtime information
	runtime = {}

	def __init__(self,path):
		#Ensure input path exisits
		if not os.path.isdir(path):
			raise ValueError("Input path: "+ path +" not found.")
		self.ids = []
		self.reads={}
		self.runtime={}
		
		for root,dirs,files in os.walk(path):
			#scan path and look for fastq files then gather ids and store in temp lists
			# print (f"fileparser path : {path}")

			for file in files:
				if '.fastq' in file or '.fastq.gz' in file:
					#get id and check if we have seen this id before by adding to id list and creating a new read object
					id = "_".join(file.split('_')[0:-1])
					
					if id not in self.ids:
						self.ids.append(id)
						self.reads[id] = self.Fastqs(id)
						print (id)
					#if fastq file is foward reads add path to .fwd
					if '_R1' in file or '_1' in file:
						if not self.reads[id].fwd:
							self.reads[id].fwd = root + '/' + file
					#if fastq file is reverese reads add path to .rev
					elif '_R2' in file or '_2' in file:
						if not self.reads[id].rev:
							self.reads[id].rev = root + '/' + file
					#if fastq file is unpaired or interleaved add path to .path
					#TODO consider the impact of this and determine best method
					else:
						if not self.reads[id].path:
							self.reads[id].paired = False
							self.reads[id].path = root + '/' + file
							
		#notify user if no fastq files were found
		if len(self.ids) == 0:
			raise ValueError("No fastq files found in " + path)

	#return a list of id with each fastq path
	def return_id_list(self):
		output_list = []
		for id in self.ids:
			if self.reads[id].paired:
				output_list.append([id,self.reads[id].fwd,self.reads[id].rev])
			else:
				output_list.append([id,self.reads[id].path])
		return output_list

	#return a list of all reads
	def return_fastq_list(self):
		output_list = []
		for id in self.ids:
			if self.reads[id].paired:
				output_list.append(self.reads[id].fwd)
				output_list.append(self.reads[id].rev)
			else:
				output_list.append(self.reads[id].path)
		return output_list

	#function to add runtime data to class
	def add_runtime(self,data_type,id,*path):
		if data_type not in self.runtime:
			self.runtime[data_type] = {}
		self.runtime[data_type][id] = path
		return None
