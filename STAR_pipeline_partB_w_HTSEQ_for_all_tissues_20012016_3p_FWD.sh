#! /bin/bash

#Karl Kremling
# April 2015
# This gets the specified SRA file from the homedir and run the tuxedo pipeline on it

#${1} is the name of the batch
#${2} is the name of the SRA file, late this should all be changed back to ${2} using ctrl + f 

batch_x=${1}
echo 'batch is' ${batch_x}
echo 'file name' ${2}
echo 'simul_processes' ${3}
echo 'cores_per_process' ${4}
echo 'total_ram' ${5}
echo 'ram_per_processor' ${6}
printf "\n"

#####FASTQDUMP SRA files to gzipped fastq
#Automate the extraction of the files from SRA format to fastq (7 minutes for 37 M reads)
#cd /workdir/${batch_x}/test_fastq_files/ 
#fastq-dump --outdir /workdir/${batch_x}/fastq_files/ --split-3 /workdir/${batch_x}/SRA_files/${2}.sra 
#rm /workdir/${batch_x}/SRA_files/${2}.sra
#
#
#
#####Randomly subsample fastq to speed up alignments
#input_dir_trim=/workdir/${batch_x}/fastq_files/
#python /workdir/${batch_x}/subsamplewithHTSeq.py 0.01 ${input_dir_trim}${2}_1.fastq ${input_dir_trim}${2}_2.fastq ${input_dir_trim}${2}_1_sub.fastq ${input_dir_trim}${2}_2_sub.fastq
#mv ${input_dir_trim}${2}_1_sub.fastq ${input_dir_trim}${2}_1.fastq 
#mv ${input_dir_trim}${2}_2_sub.fastq ${input_dir_trim}${2}_2.fastq
# 
#
#

#####TRIMMOMATIC that outputs gzipped files
#cd /workdir/${batch_x}/test_fastq_files/short_fastq
#NOTE THIS INCLUDES SE TRUseq adapter trimming
#Also includes HEADCROP 12 to remove first 12 bases
input_dir_trim=/workdir/${batch_x}/fastq_files/
java -jar /programs/trimmomatic/trimmomatic-0.32.jar SE -threads ${4} -phred33 ${input_dir_trim}${2}.fastq.gz ${input_dir_trim}trimmed/${2}_trimmed.fq ILLUMINACLIP:/programs/trimmomatic/adapters/TruSeq3-SE.fa:2:30:10 HEADCROP:12 #LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
#echo "$(du -h /workdir/${batch_x}/)" | mailx -s ${2}_job_status_workdir_trimmomatic kak268@cornell.edu

#FASTQC Qualty control on UNTRIMMED SEQUENCES 
fastqc ${input_dir_trim}${2}.fastq.gz 
mv ${input_dir_trim}${2}_fastqc /workdir/${batch_x}/fastqc/${2}_fastqc
#remove original untrimmed  fastq_files
#rm ${input_dir_trim}${2}.fastq.gz
#
#
#


#STAR alignment
outfolder=/workdir/${batch_x}/STAR_alignments_w_known_transcriptome_and_trinity_cDNAs
mkdir ${outfolder}/${2}
suboutfolder=${outfolder}/${2}
infolder=/workdir/${batch_x}/fastq_files/trimmed
/programs/STAR/STAR --genomeLoad NoSharedMemory --genomeDir /workdir/${batch_x}/STAR_index_and_fa/ --runThreadN ${4} --outReadsUnmapped Fastx --readFilesIn ${infolder}/${2}_trimmed.fq --outFilterMultimapNmax 10 --outFilterMismatchNoverLmax 0.04 --outSAMmode Full --outSAMattributes Standard --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outFileNamePrefix ${suboutfolder}/${2}_strspec_

#BWA MEM alignment
#bwa mem -t ${4} /workdir/${batch_x}/STAR_index_and_fa/ ${infolder}/${2}_trimmed.fq > ${suboutfolder}/${2}_strspec_Aligned.out.sam 

### This works: cat SRR908214_strspec_Aligned.out.sam | samtools view -bS -u -h - > SRR908214_strspec.out
# sort the output SAM file and make it a bam file do that Cufflinks can take it as input, for htseq sort by name using -n in samtools sort
#when sorting for HTseq use -n with samtools sort to sort by name
samtools view -@ ${4} -b -S -h -u -1 ${suboutfolder}/${2}_strspec_Aligned.out.sam > ${suboutfolder}/${2}_strspec_Aligned.out.bam
samtools sort -@ ${4} -m ${6}g -n ${suboutfolder}/${2}_strspec_Aligned.out.bam /workdir/${batch_x}/BAM_files/${2}_strspec.out.namesorted

#the samtools command below produces the output for cuffinks
#samtools sort -@ ${4} -m ${6}g ${suboutfolder}/${2}_strspec_Aligned.out.bam /workdir/${batch_x}/BAM_files/${2}_strspec.out.sorted

# Remove the original SAM file which was output by STAR, and the intermediate unsorted bam file
#rm ${suboutfolder}/${2}_strspec_Aligned.out.sam
rm ${suboutfolder}/${2}_strspec_Aligned.out.bam

# a way to sort using picard
#java -jar /programs/picard-tools-1.98/MergeSamFiles.jar I=${suboutfolder}/${2}_strspec_Aligned.out.sam O=/workdir/${batch_x}/BAM_files/${2}_strspec.out.sorted.bam SO=coordinate AS=false USE_THREADING=true TMP_DIR=/workdir/${batch_x}/picard_temp/temp_picard_dir_${2}/

#HTSeq quantification

python /programs/HTSeq-0.6.1/HTSeq/scripts/count.py -f bam -s yes -m union /workdir/${batch_x}/BAM_files/${2}_strspec.out.namesorted.bam /workdir/${batch_x}/known_gene_gtf/Zea_mays.AGPv3.29.gtf.gz/Zea_mays.AGPv3.29.gtf > /workdir/${batch_x}/HTSEQ_counts/${2}_HTSEQ_count.txt

#Count the mapped reads and normalize the counts by million mapped reads to get CPM
totalmappedreads=`samtools view -c -F 4 /workdir/${batch_x}/BAM_files/${2}_strspec.out.namesorted.bam` 
awk '{print $1, 1000000*$2/"'${totalmappedreads}'" }' /workdir/${batch_x}/HTSEQ_counts/${2}_HTSEQ_count.txt > /workdir/${batch_x}/HTSEQ_counts/${2}_HTSEQ_countCPM.txt 

# The quotes in the line above need to be changed for the large memory servers



#Move/rename BAM files from subfolder to main outfolder and include SRR number
#mv ${suboutfolder}/${2}_strspec.out.sorted.bam  /workdir/${batch_x}/BAM_files/${2}_strspec.out.sorted.bam 
mv ${suboutfolder}/${2}_strspec_Log.final.out ${outfolder}/${2}_strspec_Log.final.out
#rm -r ${suboutfolder}/
echo "$(du -h /workdir/${batch_x}/)" | mailx -s ${2}_job_status_workdir_STAR_alignment kak268@cornell.edu

#rm ${infolder}/${2}_1_paired.fq
#rm ${infolder}/${2}_2_paired.fq



# create a remote folder into which the cuff diff out will be placed and use iCommands and iPut to move them
#icd /iplant/home/kakremli/Maize_15DAP_Kernel/outfiles
#imkdir ${batch_x}
#icd /iplant/home/kakremli/Maize_15DAP_Kernel/outfiles/${batch_x}
#echo "#! /bin/bash" > /workdir/${batch_x}/list_of_cuffdiff_files_to_put_${batch_x}.sh  
# I would do this in a loop, but the iput command cannot take input from a loop
#ls /workdir/${batch_x}/cuffdiff_out/ | xargs -I % echo iput -r -V /workdir/${batch_x}/cuffdiff_out/${2} /iplant/home/kakremli/Maize_15DAP_Kernel/outfiles/${batch_x}/${2} >> /workdir/${batch_x}/list_of_cuffdiff_files_to_put_${batch_x}.sh 
#bash /workdir/${batch_x}/list_of_cuffdiff_files_to_put_${batch_x}.sh 

#ils | mailx -s job_status_workdir_iput kak268@cornell.edu

#icd /iplant/home/kakremli/Maize_15DAP_Kernel/outfiles
