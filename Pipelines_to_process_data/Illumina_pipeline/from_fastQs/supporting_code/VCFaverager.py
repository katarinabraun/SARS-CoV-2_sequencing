import os
import pandas as pd
from functools import reduce
from vcfClass import VCF

# This just creates a little function to move a column to the front of
# a data table, since Pandas doesn't really have a great way to do this
def movetofront(df, first):
	cols = df.columns.tolist()
	cols.remove(first)
	cols.insert(0, first)
	return (df.reindex(columns=cols))

#function to turn sample column into dataframe with each sample value given its own column, append to end of VCF dataframe
#also creates new column that contains the row's sample of origin
def extractSampleDF(df, samplename):
	Sample_column = []
	for index, row in df.iterrows():				#Goes through each row in the vcf file
		try:
			Sample_column.append(row[samplename].split(":"))	#splits up the sample column and adds it to the Sample_column list of lists
		except:
			#debug info
			print ("index = " + str(index))
			print ("row = " + str(row))
			print ("samplename = " + str(samplename))
			print ("row[samplename] = " + str(row[samplename]))

	#Turns that Sample_column list of lists into a dataframe
	Sample_column = pd.DataFrame(Sample_column, columns = df["FORMAT"][0].split(":"))
	#Adds that dataframe to the main vcf dataframe
	vcfDF = pd.concat([df.iloc[:,:-1], Sample_column], axis = 1)

	vcfDF["Sample"] = samplename			#Add a samplename column with the sample's name
	vcfDF = movetofront(vcfDF, "Sample")	#Moves that column to the front
	return vcfDF							#Finally, we return our DF

#function to reverse extractSampleDF and return to .vcf specifications
def unextractVCF(vcf, key):
	formatloc = vcf.columns.get_loc("FORMAT")
	
	vcf[key] = vcf.iloc[:,formatloc+1:].astype(str).apply(":".join, axis=1)
	print (vcf.dtypes)
	listToRemove = vcf["FORMAT"][1].split(":")
	listToRemove.append("Sample")
	for item in listToRemove.copy():
		if item not in vcf.columns:
			listToRemove.remove(item)
			
	vcf.drop(listToRemove, axis=1, inplace=True)
	return (vcf)


def VCFaverager (runCFG, repDict, vcfsToBeAveraged):
	#define default parameters
	outDir = runCFG['exec']['outdir']
	logfile = runCFG['exec']['logfile']
	outDir = os.path.join(outDir,'snp_calls')

	listofVCFs = []
	
	
	
	#Opens representative as a text file and identifies the row in which our data begins
	#vcfsToBeAveraged = [VCF(location) for location in vcfsToBeAveraged]
	print (repDict.keys())
	#go through samples, extract snp info for each replicate of that sample, merge and average freq values, then export sample-specific VCF file
	for key in repDict.keys():
		replicateVCFs = []
		for vcfname in repDict[key]:
			samplename = vcfname.split('.')[0]
			try:
				vcf = VCF(os.path.join(outDir,vcfname))
			except OSError:
				print ("\n" + os.path.join(outDir,vcfname) + " did not open appropriately. Please check file.\n")
			
			replicateVCFs.append(vcf)
			#This section splits up the "Sample_column" column from a list of 
			#terms separated by ":", assigning each one its own column.

		#merge arbitrary number of replicate DFs on SNP site (defined as same segment, same position). 
		#"reduce" applies pd.merge to two replicate DFs, then applies pd.merge to an additional DF and the newly merged DF. It repeats until only one DF remains.
		averagedVCF = reduce(lambda x, y: x.averageWithVCF(y, samplename), replicateVCFs)

		#Convert X% freq string into .X float, then averages frequency values, then rounds to 4 digits
		# mergedvcf["FREQ"] = mergedvcf["FREQ"].str.rstrip('%').astype('float')/100
		# mergedvcf["FREQ_extra"] = mergedvcf["FREQ_extra"].str.rstrip('%').astype('float')/100
		# mergedvcf["FREQ"] = round(mergedvcf[["FREQ", "FREQ_extra"]].mean(1), 4)
		
		# listOfColumnsToSum=['SDP', 'DP', 'RD', 'AD']
		# for column1 in listOfColumnsToSum:
		# 	column2 = column1+"_extra"
		# 	mergedvcf[column1] = mergedvcf[[column1, column2]].astype('int').sum(1)
		# 	print (mergedvcf[column1])
		# 	mergedvcf.drop(column2, axis=1)
		
		# mergedvcf["AD"] = round(mergedvcf["DP"]*mergedvcf["FREQ"],0).astype('int')
		
		# mergedvcf["RD"] = (mergedvcf["DP"]-mergedvcf["AD"])
		
		# #Join all the sample info into one column as specified by VCF format
		
		# mergedvcf = mergedvcf.drop(mergedvcf.filter(regex='_extra').columns, axis=1)
		
		# mergedvcf.drop(['RDF', 'RDR', 'ADF', 'ADR'], axis=1, inplace=True)

		# mergedvcf = unextractVCF(mergedvcf, key)
		
		#export vcfDF as .vcf file with appropriate metadata to top of vcf
		#easiest way I found to do this was to export a TSV with just the metadata, then append vcfDF to bottom of that same file
		# metadata = pd.Series(['##fileformat=VCFv4.2', '##source=JLL','##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\'Quality Read Depth of bases with Phred score >= 30\'>','##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\'Depth of variant-supporting bases (reads2)\'>', '##FORMAT=<ID=FREQ,Number=1,Type=String,Description=\'Variant allele frequency\'>'])

		# result = metadata.to_csv(os.path.join(outDir, key+"_averaged.vcf"), sep = '\t', index=False)
		# result = mergedvcf.to_csv(os.path.join(outDir, key+"_averaged.vcf"), mode = 'a', sep = '\t', index=False)
		averagedVCF.to_vcf(os.path.join(outDir, samplename+"_averaged.vcf"))
		listofVCFs.append(samplename+"_averaged.vcf")
		
	
	return (listofVCFs)