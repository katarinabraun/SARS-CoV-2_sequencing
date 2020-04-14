
import os, sys, pandas as pd

listofVCFs = []

def movetofront(df, first):
    cols = df.columns.tolist()
    cols.remove(first)
    cols.insert(0, first)
    return df.reindex(columns=cols)

def compareVCFs(vcf1, vcf2):
	files = [vcf1, vcf2]
	samples = []
	print (files)
	outDir = "/".join(vcf1.split("/")[:-1])
	for vcfname in files:
		
		suffix = vcfname[vcfname.find('.'):]
		if suffix != '.vcf':
			print (vcfname + ' is not a .vcf file')
			continue
		print ('----Adding ' + vcfname + '-----')
		vcf = vcfname.split("/")[-1]
		samplename = vcf.split('.')[0]
		
		rowstoskip = 0
		with open(vcfname, 'r') as TextVCF:
			for index, line in enumerate(TextVCF, 0):
				if '#CHROM' in line:
					rowstoskip = index
					break

		try:
			vcfDF = pd.read_csv(vcfname, sep='\t', skiprows=rowstoskip)
		except OSError as inst:
			print ('\n' + vcf + ' did not open appropriately. Please check file.\n')
		else:
			Sample_column = []
			for index, row in vcfDF.iterrows():
				Sample_column.append(row[samplename].split(':'))

			Sample_column = pd.DataFrame(Sample_column, columns=vcfDF['FORMAT'][0].split(':'))
			vcfDF = pd.concat([vcfDF.iloc[:, :-2], Sample_column], axis=1)
			vcfDF['Sample'] = samplename
			vcfDF = movetofront(vcfDF, 'Sample')
			
			vcfDF = vcfDF[['Sample', "#CHROM", 'POS', 'ALT', "DP", "AD", "FREQ"]]
			listofVCFs.append(vcfDF)
			samples.append(samplename)
			
	dataTable = listofVCFs[0].merge(listofVCFs[1], on=["#CHROM", "POS"], sort=True, how='outer', suffixes = ("", "_RePlow")).reset_index(drop=True)
	test = dataTable.to_csv(os.path.join(outDir, "test.csv"), index=False)
	for i, row in dataTable.iterrows():
		if row.isna().any():
			if row[0:7].isna().any():
				dataTable.iloc[i, 0] = "RePlow"
				
				dataTable.iloc[i, 3] = dataTable.iloc[i, 8]
				dataTable.iloc[i, 4] = dataTable.iloc[i, 9]
				dataTable.iloc[i, 5] = dataTable.iloc[i, 10]
				dataTable.iloc[i, 6] = dataTable.iloc[i, 11]
			else:
				dataTable.iloc[i, 0] = "Varscan"
		else:
			dataTable.iloc[i, 0] = "Both"
	dataTable = dataTable[['Sample', "#CHROM", 'POS', 'ALT', "DP", "AD", "FREQ", 'FREQ_RePlow']]

	
	banan = dataTable.to_csv(os.path.join(outDir, "comparison.csv"), index=False)
	# pivotTable = listofVCFs[0]['Sample', 'POS', 'FREQ']
	# pivotTable = dataTable.pivot(values='FREQ', index='POS', columns='Sample')
	# results = pd.ExcelWriter(outDir + '/VCFResults.xlsx')
	# dataTable.to_excel(results, 'VCF Results')
	# pivotTable.to_excel(results, 'Mutation Frequencies')
	# results.save()