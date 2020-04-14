#! /bin/bash
# exit if error occurs
set -e
## as written if we are short a bar code read count it will continue indefinatly.
## NOTE 1 of 2: All folders must has a trailing slash "/"
## NOTE 2 of 2: Edit where there is ## Manually edit this line of code 
CORES=2
## MIN_LENGTH=400
## MAX_LENGTH=700
NORMALIZE=200
## local folder where the fastq files are transfered to from grid ion
FASTQ_FOLDER="/scratch/fastq/"
## where we store intermediate files like the list of files
HOME_FOLDER="/scratch/"
## where we store demultiplexed files
BARCODE_FOLDER="/scratch/demultiplexed/"
##
TRIMMED_FOLDER="/scratch/trimmed/"
## final result
RESULT_FOLDER="/scratch/result/"
## barcodes of interest do not set higher than expected or it will not total
echo "making trimmed folder"
mkdir -p $TRIMMED_FOLDER
mkdir -p $RESULT_FOLDER
## BARCODE_COUNT must have 1 leading 0 to make it 2 digits.  (i.e. 02, 06, 10, 12, 24)
BARCODE_COUNT="16"
TARGET_READCOUNT=100000
GRID_FASTQ_PATH=###CONFIGURE 
echo "LOCAL FASTQ Folder : ${FASTQ_FOLDER}"
echo "BARCODE_FOLDER : ${BARCODE_FOLDER}"
echo "HOME_FOLDER : ${HOME_FOLDER}"
echo "TRIMMED_FOLDER : ${TRIMMED_FOLDER}"
echo "MIN_LENGTH : 400"
echo "MAX_LENGTH : 700"
echo "CORES : ${CORES}"
echo "BARCODE_COUNT : ${BARCODE_COUNT}"

while :;
do
	echo "rsync with the destination of interest"
	## Manually edit this line of code
	sshpass -p grid rsync -aP `IPADDRESSOFSEQEUNCINGMACHINE`:${GRID_FASTQ_PATH} $FASTQ_FOLDER ;
	cd ${FASTQ_FOLDER}
	ls *.fastq > ${HOME_FOLDER}fastq_files.txt
	echo "Getting fastq file list"
	## Filter run artic duplicate
	for f in `cat ${HOME_FOLDER}fastq_files.txt`; do 
		echo "testing if file list has been barcoded"
		cd ${BARCODE_FOLDER}
		if ! ls ${f::-6}* 1> /dev/null 2>&1; then
			echo "running  filter on the file to trim out short and long reads"
			/seqtk/seqtk comp ${FASTQ_FOLDER}${f} | awk '{ if (($2 >= 400) && ($2 <= 700)) { print} }' | cut --fields 1 > ${f}.list
			/seqtk/seqtk subseq ${FASTQ_FOLDER}${f} ${f}.list > ${TRIMMED_FOLDER}${f}
			# remove human reads
			
			
			cd $BARCODE_FOLDER
			echo "demultiplexing reads"
			artic demultiplex --threads $CORES ${TRIMMED_FOLDER}${f}
		fi
	done
	## now start the counting
	## Count needed to hit
	z=$TARGET_READCOUNT
	## total count to compare for the least value
	y=0
	## cd $BARCODE_FOLDER
	echo "looping through barcodes"
	for i in $(eval echo "{01..${BARCODE_COUNT}}"); do
		y=0
		echo "getting barcodee fastq file list"
		cd ${BARCODE_FOLDER}
		ls *NB${i}.fastq > ${HOME_FOLDER}multiplex_files.txt
		for bc_f in `cat ${HOME_FOLDER}multiplex_files.txt`; do 
			x=$(wc -l < "$bc_f")
			((x=x/4))
			((y=x+y))			
		done
		echo "Read count of NB${i} :  ${y}"
		if [ "$z" -gt "$y" ]; then
			z=$y
		fi
		if [ "$y" -gt $TARGET_READCOUNT ]; then
			echo "done with barcode ${i} now do minimapsteps"
			## cat barcode
			if ! ls ${RESULT_FOLDER}NB${i}.fastq 1> /dev/null 2>&1; then
				cat ${BARCODE_FOLDER}*NB${i}.fastq > ${RESULT_FOLDER}NB${i}.fastq
				artic minion --minimap2 --medaka --normalise $NORMALIZE --threads $CORES --scheme-directory /artic-ncov2019/primer_schemes --read-file ${RESULT_FOLDER}NB${i}.fastq nCoV-2019/V3 ${RESULT_FOLDER}BC${i};
				echo "workflow complete for barcode: ${i}"
			else
				echo "minimap steps already completed"
			fi
		else
			echo "continue"
		fi		
	done
	echo "Lowest Read count : ${z}"

	
	
	
	
	
	echo "sleeping for 3m"
	sleep 3m ;
done;