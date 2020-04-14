import pandas as pd
import os,sys


def importVCF(location):
	return VCF(location)

class VCF:
	def __init__ (self, location):
		self.vcfFileName = location

		with open (location, 'r') as vcffile:
			self._vcflines = vcffile.readlines()
		
		self._rowstoskip = self._getVCFStart()

		self.header = Header(self._vcflines[0:self._rowstoskip-1])

		
		self.samples = self._vcflines[self._rowstoskip].split('\t')[9:]
		self.mutations = [MutCall(row, self.samples) for row in self._vcflines[self._rowstoskip:]]
		self.SNPs = [mut for mut in self.mutations if mut.type == 'SNP']
		self.indels = [mut for mut in self.mutations if mut.type == 'insertion' or mut.type == 'deletion']
		self._hashedmuts = {mut.chrom:{int(mut1.pos):mut1 for mut1 in self.mutations if mut1.chrom == mut.chrom} for mut in self.mutations}
		
		
	def _getVCFStart(self):
		rowstoskip = 0
		#Opens stresentative as a text file and identifies the row in which our data begins
		for index, line in enumerate(self._vcflines):
			if "#CHROM" in line:
				rowstoskip = index+1
				break
		return rowstoskip

	def averageWithVCF(self,otherVCF,newSampleName=None):
		'''
		*assumes one sample per vcf*
		take intersection of mutations and average their coverages/frequencies etc.
		'''
		print (f'Averaging {self} with {otherVCF}')
		
		mutsnotinotherVCF = {(mut.chrom,mut.pos) for mut in self.mutations} #set
		
		for mut in otherVCF:
			try:
				self._hashedmuts[mut.chrom][mut.pos] = self._hashedmuts[mut.chrom][mut.pos].average(mut)
				mutsnotinotherVCF.discard((mut.chrom,mut.pos))
				if self._hashedmuts[mut.chrom][mut.pos] == None:
					mutsnotinotherVCF.add((mut.chrom,mut.pos))
			except KeyError:
				pass
			
		for chrom, pos in mutsnotinotherVCF:
			self.removemut(chrom, pos)
		
		self.renameSample(self.samples[0],newSampleName)
		#self.header = self.header.combineHeaders(otherVCF.header)
		return self

	def averageSamples(self,sampleA,sampleB,newSampleName=None):
		'''take intersection of mutations and average their coverages/frequencies etc.'''

		intersection = [mut for mut in self.mutations if mut.hasSample(sampleA) and mut.hasSample(sampleB)]
		#iterate through all mutations
		for mut in intersection:
			mut = mut.averageSamples(sampleA,sampleB)
		
		self.renameSample(sampleA, newSampleName)
		self.deleteSample(sampleB)
		return self

	def renameSample(self, origname, newName):
		'''update name of sample orignmae in both list of samples (self.samples)
		and in all mutations'''
		if newName == None:
			return 1
		for mut in self.mutations:
			mut.renameSample(origname, newName)
		self.samples = [sample if sample != origname else newName for sample in self.samples]

	def deleteSample(self, samplename):
		for mut in self.mutations:
			mut.deleteSample(samplename)
		self.samples.remove(samplename)
	
	def mergeVCFs(self,otherVCF):
		for mut in otherVCF:
			try:
				self._hashedmuts[mut.chrom][mut.pos].addSamples(mut.samples, mut)
				self.samples.extend(mut.samples)
			except KeyError:
				self.addMut(mut)
		return self

	def addMut(self,newMut):
		pass

	def removemut(self, chrom, pos):
		'''muts are stored in self.mutations and self._hashedmuts. 
		this removes all muts w/ chrom and pos from both lists.'''
		try:
			del self._hashedmuts[chrom][pos]
			self.mutations = [mut for mut in self.mutations if ((mut.chrom != chrom) or (mut.pos != pos))]
		except KeyError:
			pass

	def to_vcf(self, location):
		with open (location, 'w') as outfile:
			outfile.write(str(self.header))
			columnheader = "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	"+'\t'.join(self.samples)+'\n'
			outfile.write(columnheader)
			outfile.writelines([str(mutation)+'\n' for mutation in self.mutations])
	
	def fetchSNPs(self, chrom, start=0, end=None):
		return [mutation for mutation in self.mutations \
			if (mutation.chrom == chrom and mutation.pos > start and mutation.pos < end)]

	def __len__(self):
		return len(self.mutations)

	def __iter__(self):
		for mutation in self.mutations:
			yield mutation
	
	def __str__(self):
		return f"VCF containing {len (self.samples)} samples and {len(self.mutations)} mutation calls"


class MutCall:
	def __init__(self, row, samples):
		self._rawlist = row.split('\t')
		self.samples = samples
		self.chrom = self._rawlist[0]
		self.pos = int(self._rawlist[1])
		self.id = self._rawlist[2]
		self.ref = self._rawlist[3]
		self.alt = self._rawlist[4]
		self.qual = self._rawlist[5]
		self.filter = self._rawlist[6]
		self.info  = self._rawlist[7]
		self.format = self._rawlist[8]
		self.type = self._determinetype()
		self._sampledata = {sample:SampleMut(sample, data, self.format) for data, sample in zip(self._rawlist[9:], self.samples)}

	def _determinetype(self):
		if len(self.ref) > 1:
			return "deletion"
		elif len(self.alt) > 1:
			return "insertion"
		elif len(self.alt) == 1:
			return "SNP"
		else:
			return "translocation"
	
	def get(self, sample):
		return self._sampledata[sample]

	def averageSamples(self, sampleA, sampleB):
		self._sampledata[sampleA] = self._sampledata[sampleA].average(self._sampledata[sampleB])
		return self
	
	def renameSample(self,oldsamplename,newName):
		self._sampledata[newName] = self._sampledata.pop(oldsamplename)
		self.samples = [sample if sample != oldsamplename else newName for sample in self.samples]

	def deleteSample(self, sampletoDelete):
		self._sampledata.pop(sampletoDelete)
		self.samples.remove(sampletoDelete)

	def average(self, otherMut):
		sample = self.get(self.samples[0])
		sample2 = otherMut.get(otherMut.samples[0])
		if sample2.SDP == 0 or sample.SDP==0:
			return None
		else:
			sample = sample.average(sample2)

		tempinfo = self.info.split(';')
		for ADPloc, item in enumerate(tempinfo):
			if "ADP" in item:
				ADP = int(item.split('=')[-1])
				for otheritem in otherMut.info.split(';'):
					if "ADP" in otheritem:
						otherADP = int(otheritem.split('=')[-1])
						break
				break
		
		tempinfo[ADPloc] = str(int((ADP+otherADP)/2))
		self.info = "ADP="+";".join(tempinfo)
		return self
	
	def hasSample(self, samplename):
		return self._sampledata[samplename].exists()
	
	def addSamples(self,samplelist,newMut):
		pass
	
	def __str__(self):
		return '\t'.join([self.chrom, str(self.pos), self.id, self.ref,self.alt,self.qual,self.filter,self.info,self.format])+'\t'+'\t'.join(str(self.get(sample)) for sample in self.samples)
		
	def __iter__(self):
		for sample in self.samples:
			yield self._sampledata[sample]

class SampleMut:
	def __init__(self, samplename, data, formatinfo):
		self.name = samplename
		self.format = formatinfo
		self.other=[]
		
		for label, item in zip(self.format.split(":"), data.split(":")):
			item = item.strip('\n').strip('\"')
			
			if item == '.':
				item = 0
			elif item.isnumeric():
				item = int(item)
				
			if label == "GT":
				self.GT = item
			elif label == "GQ":
				self.GQ = item
			elif label == "SDP":
				self.SDP = item
			elif label == "DP":
				self.DP = item
			elif label == "RD":
				self.RD = item
			elif label == "AD":
				self.AD = item
			elif label == "FREQ":
				self.freqstring = item
				try:
					self.freq = round(float(item.rstrip('%'))/100,4)
				except:
					pass
			elif label == "PVAL":
				self.PVAL = item
			elif label == "RBQ":
				self.RBQ = item
			elif label == "ABQ":
				self.ABQ = item
			elif label == "RDF":
				self.RDF = item
			elif label == "RDR":
				self.RDR = item
			elif label == "ADF":
				self.ADF = item
			elif label == "ADR":
				self.ADR = item
			else:
				self.other.append((label,item))
		self._properties = [self.GT,self.GQ,self.SDP,self.DP,self.RD,self.AD,self.freqstring,self.PVAL,self.RBQ,self.ABQ,self.RDF,self.RDR,self.ADF,self.ADR]

	def average(self, otherSample):
		self.freq = round((self.freq+otherSample.freq)/2, 4)
		self.freqstring = str(self.freq*100)+"%"
		self.GQ = int(round((self.GQ+otherSample.GQ)/2, 0))
		self.SDP += otherSample.SDP
		self.DP += otherSample.DP
		oldAD = self.AD
		self.AD = int(round(self.DP*self.freq,0))
		oldRD = self.RD
		self.RD = self.DP-self.AD
		self.RBQ = int(round((self.RBQ+otherSample.RBQ)/2, 0))
		self.ABQ = int(round((self.ABQ+otherSample.ABQ)/2, 0))

		try:
			RDFtoRD = ((self.RDF/oldRD)+(otherSample.RDF/otherSample.RD))/2
			self.RDF = int(round(self.RD*RDFtoRD,0))
			self.RDR = self.RD-self.RDF
		except ZeroDivisionError:
			self.RDF=0
			self.RDR=0

		try:
			ADFtoRD = ((self.ADF/oldAD)+(otherSample.ADF/otherSample.AD))/2
			self.ADF = int(round(self.AD*ADFtoRD,0))
			self.ADR = self.AD-self.ADF
		except ZeroDivisionError:
			self.RDF=0
			self.RDR=0
		

		self._properties = [self.GT,self.GQ,self.SDP,self.DP,self.RD,self.AD,self.freqstring,self.PVAL,self.RBQ,self.ABQ,self.RDF,self.RDR,self.ADF,self.ADR]
		
		print (self._properties)
		return self

	def exists(self):
		if self.SDP == 0:
			return False
		else:
			return True
	
	def __str__(self):
		return ":".join([str(item) for item in self._properties])

class Header:
	def __init__(self, headertext):
		self.text = headertext

	def combineHeaders(self,otherheader):
		for line in otherheader.text:
			if line not in self.text:
				self.text.append(line)
	
	def __str__(self):
		return "".join(self.text)