import os
import calldocker as cd
from shutil import copyfile

#Utilizes bcftools to index, compress, and combine a list of vcf files using bcftools
def VCFcombiner(runCFG, vcffiles, outputname):
	
	if vcffiles[0].split('/')[:-2] != runCFG['exec']['outdir']:
		workingDir = '/'.join(vcffiles[0].split('/')[:-1])
	elif runCFG['exec']['outdir'].split("/")[-1] != 'snp_calls':
		workingDir = os.path.join(runCFG['exec']['outdir'],'snp_calls')
	else:
		workingDir = os.path.abspath(runCFG['exec']['outdir'])
	if len(vcffiles)<2:
		copyfile (os.path.join(workingDir, vcffiles[0]), os.path.join(workingDir, "all_snps.vcf"))
		return ("all_snps.vcf")
	cmds = []
	listofVCFs = []
	if len(vcffiles) > 1:
		#bcftools requires all vcfs to be indexed and compressed. This for loop generates commands to do that for every vcf in vcffiles.
		for vcf in vcffiles:
			vcflocal = vcf.split("/")[-1]
			cmd = (f"bgzip -i {vcflocal} && tabix -p vcf {vcflocal}.gz")
			cmds.append(cmd)
			listofVCFs.append(vcflocal+".gz")

		#generate command for bcftools merge
		master_vcf_input_list = " ".join(listofVCFs)
		cmds.append(f"bcftools merge -m none {master_vcf_input_list} > {outputname}")
		cmd = " && ".join(cmds)
		cmds = ['sh', '-c', cmd] #-sh and -c allow && to work in a docker container
		cd.call(cmds, '/snpcalls', {workingDir:"/snpcalls"})

	elif len(vcffiles) < 1:
		print ("No VCFs to merge!\n\n")
	#if asked to merge one vcf file, just return that vcf file
	elif len(vcffiles) == 1:
		outputname = vcffiles[0].split("/")[-1]
	return(outputname)