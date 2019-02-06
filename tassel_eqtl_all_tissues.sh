#!/bin/bash

#---------------------------------------------------------------------
# Title........... Karl TASSEL script for Nature 2018 results
# Author.......... Brandon Monier
# Created......... 2018-10-08 16:16:12 EST
# Last Modified... 2018-10-10 15:57:16 EST
#
#
# Objective:
#   - This script will run Karl's Nature 2018 data quickly and
#     iterate over all 10 maize chromosomes
#
#   - This script will also iterate through each of Karl's tissue
#     traits
#---------------------------------------------------------------------

# Assign Karl's tissues to variable for iteration
traits_all="GShoot Kern L3Base L3Tip LMAD LMAN"

# Iterate through it all...
for trait in $traits_all
do
    ## Copy expression matrices to working directory
    printf "\n#------\n# Copying ${trait} expression matrix to /workdir...\n#------\n\n"
    cp /home/bm646/tassel_runs/karl_expression_matrix/TASSEL_HEADER_df_STAR_HTSeq_counts_B73_match_based_on_genet_dist_DESeq2_normed_rounded_origNames_and_Altnames_${trait}_sm_rand_and_box_coxed_w_MDSPCs_and_25PEERs_avgDups.txt /workdir/bm646/

    for chrom in {1..10}
    do
        ### Copy files to working directory
        echo "Copying chromosome ${chrom} to /workdir..."
        cp /home/bm646/tassel_runs/karl_chrom_hmp/merged_flt_c${chrom}.hmp321.onlyRNAset_MAFover005.KNNi.hmp.txt /workdir/bm646/

        ### TASSEL Run
        echo "Running TASSEL 5 on chromosome ${chrom} for trait ${trait}..."
        /programs/tassel-5-standalone/run_pipeline.pl \
            -debug /workdir/bm646/debug_logs/${trait}_chr${chrom}.log \
            -Xmx1000g \
            -maxThreads 100 \
            -fork1 \
            -r /workdir/bm646/TASSEL_HEADER_df_STAR_HTSeq_counts_B73_match_based_on_genet_dist_DESeq2_normed_rounded_origNames_and_Altnames_${trait}_sm_rand_and_box_coxed_w_MDSPCs_and_25PEERs_avgDups.txt \
            -fork2 -h /workdir/bm646/merged_flt_c${chrom}.hmp321.onlyRNAset_MAFover005.KNNi.hmp.txt \
            -combine3 \
            -input1 \
            -input2 \
            -intersect \
            -FastMultithreadedAssociationPlugin \
            -MaxPValue 0.00001 \
            -writeToFile true \
            -outputFile /workdir/bm646/${trait}_output/${trait}_merged_flt_c${chrom}_maxP000001.hmp321.onlyRNAset_MAFover005.KNNi.w.df_STAR_HTSeq_counts_B73_match_based_on_genet_dist_DESeq2_normed_rounded_origNames_and_Altnames_sm_rand_and_box_coxed_w_MDSPCs_and_25PEERs_avgDups.txt \
            -runfork1 \
            -runfork2 \

        echo "Finished analyzing chromosome ${chrom} for trait ${trait}..."
    done
done

echo "Quitting program..."
/programs/bin/labutils/endres.pl
