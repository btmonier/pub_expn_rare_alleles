#Karl Kremling
#October 2014

#Practicing with PEER normalization

library(peer)
#install.packages("readr")
library(readr)
#library(qtl)
#cross <- read.cross(format="csvs", genfile="SNP_single_15DAP_SNP_for_PEER_tpose.csv", phefile="CPMs_log_w_rand_smalladdition_from_STAR_alignmentof_all_368_15DAP_RNAseq_15DAP_no_txt_in_header_tpose_header_now_ID.csv", genotypes=c(0,2))

covariates=TRUE
PEER_number=25
#i="df_STAR_HTSeq_counts_B73_match_based_on_genet_dist_DESeq2_normed_rounded_origNames_and_Altnames_sm_rand_and_box_coxed_L3Tip.txt"
#i="df_STAR_HTSeq_counts_B73_match_based_on_genet_dist_DESeq2_normed_rounded_origNames_and_Altnames_noL3Mid_sm_rand_and_box_coxed.txt"
basedir="/media/kak268/A_2TB_Internal/RNA/Expressions_quants/HTSEQ_counts_STAR_against_Zea_mays_B73_AGPv3.29_w_gtf_annotation/count_mat/count_mats_matched_to_AltNames/count_mats_separated_by_tissue_and_individually_box_coxed/"
#for (i in list.files(basedir, pattern=".txt")){
#expression_df=read_delim(file=paste(basedir, i, sep=""), col_names=T, delim=" ")
for (i in c("GRoot","GShoot","Kern","L3Base","L3Tip","LMAD","LMAN")){ #allnonL3Mid, for all taxa
cat(i, "\n")
  
expression_df=read_delim(file=paste(basedir, "df_STAR_HTSeq_counts_B73_match_based_on_genet_dist_DESeq2_normed_rounded_origNames_and_Altnames_sm_rand_and_box_coxed_", i,".txt", sep=""), col_names=T, delim=" ")

Taxa_list_w_tissue=expression_df[,1:10]
colnames(expression_df)[6]="ID"
expression_df_only_ID=expression_df[-c(1:5, 7:10)]
Taxa_list=as.data.frame(expression_df_only_ID$ID)
colnames(Taxa_list)="HMP32Name"

#Accoutn for the Genetic PCs when calculating the PEERs
if (covariates==TRUE){
covdir="/media/kak268/B_2TB_Internal/Genotypes/HM32/LDKNNi/PCs/"
#covs=read.csv(paste(covdir, "PCs_for_PEER_factor_calculation_merged_flt_c1-10.hmp321.onlyL3TipRNAset_2pctsubset.KNNi.hmp.csv", sep=""), header=T) # make sure these covariates are the same order as the 
covs=read.table(paste(covdir, "MDS_5coords_PEER_fmt_merged_flt_c1-10.hmp321.onlyRNAset_MAFover005.KNNi.txt", sep=""), header=T )
#repeat rows from the covariate table to match the phenotype length so that the taxa are repeated which have multiple expression phenotypes or tissues
covs_w_rep=merge(x=Taxa_list, y=covs, by="HMP32Name")
}


model=PEER()

PEER_setPhenoMean(model, as.matrix(expression_df_only_ID[,2:ncol(expression_df_only_ID)]))
dim(PEER_getPhenoMean(model))

PEER_setNk(model, PEER_number)

if (covariates==TRUE){
  PEER_setCovariates(model, as.matrix(covs_w_rep[, 2:ncol(covs_w_rep)]))
  dim(PEER_getCovariates(model)) 
}

num_of_PCs=ncol(PEER_getCovariates(model))

PEER_update(model) # if you accidentally import a different number of cov's and try to run PEER then R will crash 
PEER_plotModel(model)


#Write the PEER factors to a file
#tposed_covs_w_PEER_factors=t(cbind(as.matrix(covs_w_rep[,1]),(as.matrix(PEER_getX(model)))))
covs_w_PEER_factors=as.matrix(PEER_getX(model))

cov_name_list="HMP32Name"
for (x in 1:ncol(covs_w_PEER_factors)){ #add row names of cov or PEER factor
  if (x<=num_of_PCs){
    cov_name_list=c(cov_name_list, paste("geneticPC",x,sep=""))
  }
  else # subtract num_of_PCs from x below becasue there are 5 covariates ahead of the peer factors
    cov_name_list=c(cov_name_list, paste("PEER_factor",x-num_of_PCs,sep=""))
}

covs_w_PEER_factors_and_TaxaNames=cbind(Taxa_list, covs_w_PEER_factors)
colnames(covs_w_PEER_factors_and_TaxaNames)=cov_name_list

covs_w_PEER_factors_and_TaxaNames_w_TissueNames=cbind(Taxa_list_w_tissue, covs_w_PEER_factors_and_TaxaNames[,2:ncol(covs_w_PEER_factors_and_TaxaNames)])

outputdir="/media/kak268/B_2TB_Internal/Genotypes/HM32/LDKNNi/PCs/"
write.table(covs_w_PEER_factors_and_TaxaNames_w_TissueNames, file=paste(outputdir, i,"/", "MDSPCs_from_HMP321_w_", PEER_number,"_PEER_factors_", i,"_complete_taxa_and_tiss_names.txt", sep=""), quote=F, sep="\t", row.names=F, col.names=T)



#Write a version of the covs, PEERs and the expression as one DF to be used in TASSEL
covs_w_PEER_factors_and_TaxaNames_w_Expression_for_TASSEL=cbind(Taxa_list, covs_w_PEER_factors_and_TaxaNames[,2:ncol(covs_w_PEER_factors_and_TaxaNames)], expression_df_only_ID[,2:ncol(expression_df_only_ID)])

colnames(covs_w_PEER_factors_and_TaxaNames_w_Expression_for_TASSEL)[1]="Taxa"
TASSEL_filename_expr_w_covs_PEERs=paste(basedir, "expression_df_w_covs/", i,"/", "df_STAR_HTSeq_counts_B73_match_based_on_genet_dist_DESeq2_normed_rounded_origNames_and_Altnames_", i ,"_sm_rand_and_box_coxed_w_MDSPCs_and_", PEER_number,"PEERs.txt", sep="")
write.table(covs_w_PEER_factors_and_TaxaNames_w_Expression_for_TASSEL, file=TASSEL_filename_expr_w_covs_PEERs, quote=F, sep="\t", row.names=F, col.names=T)



TASSEL_preheader_for_covs_and_exp_values=c("taxa", rep("covariate", ncol(covs_w_PEER_factors)), rep("data", ncol(expression_df_only_ID)-1))
# after writing this to fil, manually add <Phenotype> to the top of the preheader file and then use cat to paste it above the expression values
preheader_filename=paste(basedir, "expression_df_w_covs/", i,"/", i,"_TASSEL_preheader_for_covs_and_exp_values", sep="")
write.table(t(TASSEL_preheader_for_covs_and_exp_values), file=preheader_filename, quote=F, sep="\t", row.names=F, col.names=F)
system(paste("sed '1i <Phenotype>' ", preheader_filename, TASSEL_filename_expr_w_covs_PEERs, ">", paste(basedir, "expression_df_w_covs/", i,"/", "TASSEL_HEADER_df_STAR_HTSeq_counts_B73_match_based_on_genet_dist_DESeq2_normed_rounded_origNames_and_Altnames_", i ,"_sm_rand_and_box_coxed_w_MDSPCs_and_", PEER_number,"PEERs.txt", sep=""), sep=" "))
}
#system(paste("sed '1i <Phenotype>' ", preheader_filename, TASSEL_filename_expr_w_covs_PEERs, ">", paste(basedir, "expression_df_w_covs/", i,"/", "TASSEL_HEADER_df_STAR_HTSeq_counts_B73_match_based_on_genet_dist_DESeq2_normed_rounded_origNames_and_Altnames_", i ,"_sm_rand_and_box_coxed_w_MDSPCs_and_", PEER_number,"PEERs.txt", sep=""), sep=" "))