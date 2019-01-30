# README


**This repository holds the scripts which are associated with the publication "Dysregulation of expression corelates with rare allele burden and fitness loss in maize".**
**Sequence data that support the findings of this study have been deposited in the SRA under SRP115041/PRJNA383416 at https://www.ncbi.nlm.nih.gov/bioproject/383416.**
**Processed expression counts are available at the Cyverse Discovery Environment (de.cyverse.org/de/) under directory: /iplant/home/shared/panzea/dataFromPubs/.**


**Processing of RNASeq reads**
=============================
#### Trimming (Trimmomatic), alignment (STAR), and counting(HTSeq) partA (This sets up the indices for alignment and parallelizes over alignments done in part B)
STAR_pipeline_partA_w_HTSEQ_for_all_tissues_20012016_3p_FWD.sh
#### Trimming (Trimmomatic), alignment (STAR), and counting(HTSeq) partB
STAR_pipeline_partB_w_HTSEQ_for_all_tissues_20012016_3p_FWD.sh
#### Post alignment processing of expression counts with DESEQ2
DESeq_normalization_of_HTSEQ_counts.R



**Figure 1, Extended Data Figure 2, 3, 4**
======
#### Average duplicates and format DESeq normalized counts (on CyVerse) and make tissue-wise expression dataframe for use in downstream scripts
Extract_genes_from_each_tissue_AvgOverDups.R
#### Make numerical rare allele hapmap file that contains variants only near genes
Extract_gene_proximal_rare_alleles_and_numericalize.R
#### Script that extracts Feature_Name, Taxon, and Upstr_5kb_rare_homo_ct from numericalized genotypes of near gene rare alleles
Extract_GeneName_TaxaName_RareCt_for_use_in_Each_tissue_vs_rare_alleles_5000genes_and_all_genes_v5.R
#### Correlate expression rank with local rare allele abundance and make smile plots
Each_tissue_vs_rare_alleles_5000genes_and_all_genes_v5_for_paper.R


**Figure 2, Extended Data Figure 5, 6**
======
#### Run TASSEL Version: 5.2.32 fast association plugin for eQTL mapping
multi_thread_assoc_agpv3.sh
#### Filter eQTL output from TASSEL for hits in a cis window
Filter_eQTL_from_TASSEL_for_cis.R
#### Pull out the top hits in a certain window from eQTL results
Record_top_hit.R
#### Plot 2D MAF plots colored by eQTL R2
MAF_comparison_between_Tropical_and_RNAset_and_eQTL_hits_CBSU_v5_cis_only.R

**Figure 3, Extended Data Figure 7**
======
#### Extract top 5000 genes and avg over duplicates
Extract_genes_from_each_tissue_AvgOverDups.R
#### Ridge regression to predict phenotype from expression and expression dysregulation
expr_deviation_vs_fitness_5000genes_all_tissue_ridgeregression_v4_nested_cross_validation.R


**Ext Data Figure 8, Ext data table 1**
====================
#### Correlate dysregulation with fitness
expr_deviation_vs_fitness_5000genes_all_tissue_v5.R


**Hidden Factor calculation**
===============
#### PEER Hidden Factors
PEER_calculation.R

