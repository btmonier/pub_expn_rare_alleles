#! /bin/bash

#Karl Kremling
# April 2015
# This gets the SRA files from the homedir at CBSU and runs the STAR + Cufflinks pipeline in parallel on 10 files at a time


# use: STAR_pipeline_partA.sh

cd /workdir/
#Decide name of the directory, and command file b and num of cores for this batch
submission_name=all_tissues_20012016_3p_FWD

batch_x=kak268/${submission_name}
command_file_a=STAR_pipeline_partA_w_HTSEQ_for_${submission_name}.sh
command_file_b=STAR_pipeline_partB_w_HTSEQ_for_${submission_name}.sh
#cores=`grep -c ^processor /proc/cpuinfo`
num_of_processors=128 #`echo ''${cores}'/'2'' | bc`  #DIVIDE by 2 if processors are dual core
total_ram=1000 #$(free -g | awk '/^Mem:/{print $2}')

cores_per_process=4

ram_per_processor=`echo ''${total_ram}'/'${num_of_processors}'' | bc`
simul_processes=`echo ''${num_of_processors}'/'${cores_per_process}'' | bc`

echo 'user specified values:'
printf "\n"
echo 'batch_name =' ${batch_x}
echo 'command_file_b =' ${command_file_b}
echo 'num_of_processors =' ${num_of_processors}
echo 'total_ram =' ${total_ram}
echo 'ram_per_processor =' ${ram_per_processor}
echo 'simul_processes =' ${simul_processes}
echo 'cores_per_process =' ${cores_per_process}
printf "\n\n\n"

#Make the directories which will contain the intermediate files
#mkdir /workdir/kak268
mkdir /workdir/${batch_x}
#mkdir /workdir/${batch_x}/SRA_files
mkdir /workdir/${batch_x}/parallelize_folder
mkdir /workdir/${batch_x}/fastq_files
#mkdir /workdir/${batch_x}/fastq_files/trimmed_paired
#mkdir /workdir/${batch_x}/fastq_files/trimmed_unpaired
mkdir /workdir/${batch_x}/fastq_files/trimmed
mkdir /workdir/${batch_x}/STAR_index_and_fa
mkdir /workdir/${batch_x}/BWAMEM_index_and_fa
mkdir /workdir/${batch_x}/known_gene_gtf
mkdir /workdir/${batch_x}/STAR_known_gene_index
mkdir /workdir/${batch_x}/STAR_alignments_w_known_transcriptome_and_trinity_cDNAs
mkdir /workdir/${batch_x}/picard_temp/
mkdir /workdir/${batch_x}/BAM_files
mkdir /workdir/${batch_x}/HTSEQ_counts
mkdir /workdir/${batch_x}/fastqc
#mkdir /workdir/${batch_x}/cuffdiff_out/

#Copy the perl parallelization script 
cp /programs/bin/perlscripts/perl_fork_univ.pl /workdir/${batch_x}/parallelize_folder/perl_fork_univ.pl


###Copy the command files in to the analysis directory
#note to self: the fastq download file can be changed so file names make sense by using this command which changes the names of the fastq files to reflect the taxon and tissue name: 
# sed 's/5664.*KAKRNA11_//' download_miseq_10343057.sh > download_10343057_newnames.sh

cp /workdir/kak268/${command_file_a} /workdir/${batch_x}/${command_file_a}
cp /workdir/kak268/${command_file_b} /workdir/${batch_x}/${command_file_b}
#cp /workdir/kak268/download_${submission_name}_newnames.sh /workdir/${batch_x}/download_${submission_name}_newnames.sh

#Copy all the fastq files from the core to the local disk
cd /workdir/${batch_x}/fastq_files
###THIS IS NORMALLY WHERE THE WGET COMMANDS FOR THE FASTQ files would go, or where the bash download script is used
# Bash command to download files
#perl /workdir/${batch_x}/parallelize_folder/perl_fork_univ.pl /workdir/${batch_x}/download_${submission_name}_newnames.sh 6
#cp command to copy fastq files into the batch directory
cp /workdir/kak268/fastq_repository/*.fastq.gz /workdir/${batch_x}/fastq_files


cd /workdir/

##### Copy maize refgen files needed for pipeline
#COPY index from CBSU shared folder
#remove symbolic link to maize.fa then Copy whole genome fasta file maize.fa
# copy the refgen fasta file from the shared folder (no longer needed since it is inside of the STAR index folder which I copy from my drive
#Copy .gtf annotations from ensemble
#wget ftp://ftp.ensemblgenomes.org/pub/plants/release-23/gtf/zea_mays/Zea_mays.AGPv3.23.gtf.gz -P /workdir/${batch_x}/known_gene_gtf/Zea_mays.AGPv3.23.gtf.gz
# may need to use the -f flag in the gunzip to overwrite teh zipper file
#gunzip /workdir/${batch_x}/known_gene_gtf/Zea_mays.AGPv3.23.gtf.gz/Zea_mays.AGPv3.23.gtf.gz # > /workdir/${batch_x}/known_gene_gtf/Zea_mays.AGPv3.23.gtf.gz/Zea_mays.AGPv3.23.gtf
#
#####Create the STAR known transcriptome-index # didn't actually use trinity stuff, but directory name is still used
#/programs/STAR_2.4.2a/STAR --runMode genomeGenerate -runThreadN 24 --genomeDir transcriptome/STAR_cDNA_idx/ --genomeFastaFiles transcriptome/Zea_mays.AGPv3.29.cdna.all_w_Trinity_exc_NonAngiospermae_exc_over_85pct_B73_hits.fa --limitGenomeGenerateRAM 120000000000

#####Create the bwa mem known transcriptome index
#bwa index ../../transcriptome/Zea_mays.AGPv3.29.cdna.all_w_Trinity_exc_NonAngiospermae_exc_over_85pct_B73_hits.fa

#Copy Star index from AGPv3.29 + Trinity assembly excluding non angiospermidae hits and >85% B73 cDNA hits
#cp -r /home/kak268/Desktop/kak268/STAR_genome_indices/STAR_idx_Zea_mays.AGPv3.29.cdna.all_w_Trinity_exc_NonAngiospermae_exc_over_85pct_B73_hits/* /workdir/${batch_x}/STAR_index_and_fa

#Copy BWAMEM index from AGPv3.29 + Trinity assembly excluding non angiospermidae hits and >85% B73 cDNA hits
#cp -r /home/kak268/Desktop/kak268/BWAMEM_genome_indices/BWAMEM_idx_Zea_mays.AGPv3.29.cdna.all_w_Trinity_exc_NonAngiospermae_exc_over_85pct_B73_hits/* /workdir/${batch_x}/BWAMEM_index_and_fa



##### Copy maize refgen files needed for pipeline
#COPY index from CBSU shared folder
#remove symbolic link to maize.fa then Copy whole genome fasta file maize.fa
# copy the refgen fasta file from the shared folder (no longer needed since it is inside of the STAR index folder which I copy from my drive
#Copy .gtf annotations from ensemble
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-29/gtf/zea_mays/Zea_mays.AGPv3.29.gtf.gz -P /workdir/${batch_x}/known_gene_gtf/Zea_mays.AGPv3.29.gtf.gz
# may need to use the -f flag in the gunzip to overwrite teh zipper file
gunzip /workdir/${batch_x}/known_gene_gtf/Zea_mays.AGPv3.29.gtf.gz/Zea_mays.AGPv3.29.gtf.gz 
echo "$(du -h /workdir/${batch_x}/)" | mailx -s job_status_workdir_Fasta_and_STARindex_copied kak268@cornell.edu

#
#
#
#####Create the STAR known transcriptome-index
#/programs/STAR/STAR --runMode genomeGenerate --runThreadN 16 --genomeDir /workdir/${batch_x}/STAR_index_and_fa/ --genomeFastaFiles /workdir/${batch_x}/STAR_index_and_fa/maize.fa --sjdbOverhang 89 --sjdbGTFfile /workdir/${batch_x}/known_gene_gtf/Zea_mays.AGPv3.23.gtf.gz/Zea_mays.AGPv3.23.gtf
#/programs/STAR/STAR --runMode genomeGenerate --runThreadN 64 --genomeDir /workdir/kak268/transcriptome --genomeFastaFiles /workdir/${batch_x}/STAR_index_and_fa/maize.fa --sjdbOverhang 89 --sjdbGTFfile Zea_mays.AGPv3.29.gtf.gz/Zea_mays.AGPv3.29.gtf

cp /home/kak268/Desktop/kak268/STAR_genome_indices/maizerefgenv3_Zea_maysAGPv3_29_gtf/* /workdir/${batch_x}/STAR_index_and_fa/
#cp /workdir/batch_test6_old/STAR_index_and_fa/* /workdir/${batch_x}/STAR_index_and_fa/
#


echo "$(du -h /workdir/${batch_x}/)" | mailx -s job_status_workdir_Fasta_and_STARindex_copied kak268@cornell.edu

#
#


echo "$(du -h /workdir/${batch_x}/)" | mailx -s job_status_workdir_STAR_index_and_fa kak268@cornell.edu
#
#
#
#####Create list of all fastq/SRA files to feed to _pipeline_partB.sh
#ls /home/kak268/Desktop/kak268/kernel_Fu_SRA_files | grep 'SRR908' | sed 's/_.*//' | xargs -I % echo bash /workdir/${command_file_b} ${batch_x} ${1} % ${simul_processes} ${cores_per_process} ${total_ram} ${ram_per_processor} > calls_to_${command_file_b}.txt
ls /workdir/${batch_x}/fastq_files | grep '.fastq.gz' | sed 's/.fastq.gz//'| xargs -I % echo bash /workdir/${batch_x}/${command_file_b} ${batch_x} % ${simul_processes} ${cores_per_process} ${total_ram} ${ram_per_processor} > /workdir/${batch_x}/calls_to_${command_file_b}.txt

perl /workdir/${batch_x}/parallelize_folder/perl_fork_univ.pl /workdir/${batch_x}/calls_to_${command_file_b}.txt ${simul_processes}

grep 'Total Sequences' /workdir/${batch_x}/fastqc/*/fastqc_data.txt | perl -pe 's/\/workdir.*?(fastqc\/)//g' | sed 's/_fastqc.*:/ /' | sed 's/Total Sequences/Total_Sequences/' > /workdir/${batch_x}/${submission_name}count_of_input_sequences.txt

echo "$(du -h /workdir/${batch_x}/)" | mailx -s job_status_workdir_all_alignments_done kak268@cornell.edu

#


#####Move all BAM files to new folder 
### this is now done automatically inside of the ...partB
#cp /workdir/${batch_x}/STAR_alignments_w_known_transcriptome_and_trinity_cDNAs/SRR*/*.out.sorted.bam /workdir/${batch_x}/BAM_files/
#
#
#
#
#####cuffdiff expression
#ls /workdir/${batch_x}/BAM_files/*.sorted.bam  > /workdir/${batch_x}/bam_list.txt
#cuffdiff -o /workdir/${batch_x}/cuffdiff_out/ -p ${num_of_processors} --total-hits-norm --library-norm-method classic-fpkm -b /workdir/${batch_x}/STAR_index_and_fa/maize.fa -u /workdir/${batch_x}/known_gene_gtf/Zea_mays.AGPv3.23.gtf.gz/Zea_mays.AGPv3.23.gtf -library-type fr-firststrand $(cat /workdir/${batch_x}/bam_list.txt)
#echo "$(du -h /workdir/${batch_x}/)" | mailx -s job_status_workdir_cuffdiff_all_complete kak268@cornell.edu
#
#
#
#
#
#
