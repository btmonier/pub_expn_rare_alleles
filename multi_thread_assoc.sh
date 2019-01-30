#!/bin/sh
#Work with avgDups Files
# this script is passed to a loop and calls FastMultithreadedAssociationPlugin
# this script then loops over the chromosomes one by one
#cat tissue_list.txt | while read line ; do bash multi_thread_assoc.sh $line ; done

#copy the following into tissue_list.txt
#LMAD
#LMAN
#GRoot
#GShoot
#L3Tip
#L3Base
#Kern

#if not already done, you need to download tassl5 via git into your local /workdir
#git clone https://bitbucket.org/tasseladmin/tassel-5-standalone.git

#export JAVA 8 the PATH so that tassel5 can run
export JAVA_HOME=/usr/local/jdk1.8.0_45;
export PATH=$JAVA_HOME/bin:$PATH;
tissue=$1
echo starting ${tissue} tissue

#perl /workdir/kak268/tassel-5-standalone/run_pipeline.pl -debug /workdir/kak268/debug_logs/${tissue}.log -Xmx500g -maxThreads 63 -fork1 -r /workdir/kak268/phenotypes/TASSEL_HEADER_df_STAR_HTSeq_counts_B73_match_based_on_genet_dist_DESeq2_normed_rounded_origNames_and_Altnames_${tissue}_sm_rand_and_box_coxed_w_MDSPCs_and_25PEERs.txt -fork2 -h /workdir/kak268/merged_flt_c1-10.hmp321.onlyRNAset_MAFover005.KNNi.hmp.txt -combine3 -input1 -input2 -intersect -FastMultithreadedAssociationPlugin -MaxPValue 0.000001 -writeToFile true -outputFile /workdir/kak268/results/${tissue}_merged_flt_c1-10.hmp321.onlyRNAset_MAFover005.KNNi.w.df_STAR_HTSeq_counts_B73_match_based_on_genet_dist_DESeq2_normed_rounded_origNames_and_Altnames_sm_rand_and_box_coxed_w_MDSPCs_and_25PEERs.txt -runfork1 -runfork2

#echo finished ${tissue} tissue
#echo "$(du -h /workdir/kak268/)" | mailx -s eQTL_GWAS_finished_${tissue}_tissue kak268@cornell.edu


for i in {1..10}; do 
perl /workdir/kak268/tassel-5-standalone/run_pipeline.pl -debug /workdir/kak268/debug_logs/${tissue}chr${i}.log -Xmx500g -maxThreads 92 -fork1 -r /workdir/kak268/phenotypes/TASSEL_HEADER_df_STAR_HTSeq_counts_B73_match_based_on_genet_dist_DESeq2_normed_rounded_origNames_and_Altnames_${tissue}_sm_rand_and_box_coxed_w_MDSPCs_and_25PEERs_avgDups.txt -fork2 -h /workdir/kak268/genotypes/merged_flt_c${i}.hmp321.onlyRNAset_MAFover005.KNNi.hmp.txt -combine3 -input1 -input2 -intersect -FastMultithreadedAssociationPlugin -MaxPValue 0.00001 -writeToFile true -outputFile /workdir/kak268/results/${tissue}_merged_flt_c${i}_maxP000001.hmp321.onlyRNAset_MAFover005.KNNi.w.df_STAR_HTSeq_counts_B73_match_based_on_genet_dist_DESeq2_normed_rounded_origNames_and_Altnames_sm_rand_and_box_coxed_w_MDSPCs_and_25PEERs_avgDups.txt -runfork1 -runfork2

echo finished ${tissue} tissue chr $i
echo "$(ls -lht /workdir/kak268/results)" | mailx -s eQTL_GWAS_finished_${tissue}_tissue_chr_${i} kak268@cornell.edu ; done

