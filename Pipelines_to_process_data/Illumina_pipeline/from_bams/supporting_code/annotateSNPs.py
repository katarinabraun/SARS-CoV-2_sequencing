import os
import shutil

def annotateSnpeff(runCFG, ref, vcf, output_dir, sample_name, threads):

	output_dir = os.path.join(runCFG['exec']['outdir'], 'snp_calls')
	checkexists(output_dir)

	# run snpEff to annotate vcf file
	# only annotate variants within features (by setting ud = 0)
	print('--Annotate variant impact on coding sequence with snpEff--')

	annotated_vcf = output_dir + "/" + sample_name + "_annotated.vcf"

	with open(annotated_vcf, "wb") as out:
		snpeff_cmd = ['java',
					  '-jar',
					  appRoot + '/bin/snpEff.jar',
					  '-ud',
					  '-onlyProtein'
					  '0',
					  ref,
					  vcf]

		subprocess.call(snpeff_cmd, stdout=out)

	# replace generic Sample1 sample identifier with sample_name
	with fileinput.FileInput(annotated_vcf, inplace=True) as file:
		for line in file:
			print(line.replace('Sample1', sample_name), end='')

	# copy snpEff summary file to output location
	shutil.copyfile('snpEff_summary.html', output_dir + '/log/snpEff_summary.html')

	# remove snpEff summary file and text summary invokation location
	os.remove('snpEff_summary.html')
	os.remove('snpEff_genes.txt')

	return snpeff_cmd