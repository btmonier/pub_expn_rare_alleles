#Karl Kremling
#Aoril 3 2016
#extract top 5000 genes for each tissue withaveraging over duplicate taxa
#Separate out single tissue of untrnsofrmed, but DESeq normalized, counts NOT BOX COXED SINCE THAT WOULD CAUSE LOSS OF ability to compare between genes

#install.packages("miscTools")
library("miscTools")
library("readr")
library("reshape")
library("plyr")
library("doParallel")
#install.packages("registerDoMc")
#library("registerDoMc")
#install.packages("doMC")
library("doMC")
registerDoMC(cores=12)


#
#
#
basedir="/media/kak268/A_2TB_Internal/RNA/Expressions_quants/HTSEQ_counts_STAR_against_Zea_mays_B73_AGPv3.29_w_gtf_annotation/count_mat/count_mats_matched_to_AltNames/"
expression_df=read.table(paste(basedir, "df_STAR_HTSeq_counts_B73_match_based_on_genet_dist_DESeq2_normed_rounded_origNames_and_Altnames.txt", sep=""), header=T)


gff_file=read.delim("/home/kak268/work_data/maize_genome/Zea_mays_gene_only.AGPv3.29.gff3", header=F, sep="\t", quote="") # add skip=6 if it's a full gff file
colnames(gff_file)=c("Chr", "Annot_source", "Feature_type", "Start", "End", "blankcol1", "Strand", "blankcol2", "Feature_Name")
transcript_file=gff_file[gff_file$Feature_type=='gene',] # keep only the entries in the GFF file which are transcripts/genes
transcript_file$Feature_Name=sub(pattern= "ID=gene:", replacement='', x=transcript_file$Feature_Name)
transcript_file$Feature_Name=sub(pattern= ";.*", replacement='', x=transcript_file$Feature_Name)
transcript_file_only_chr1_to_chr10_gene_list=transcript_file[transcript_file$Chr %in% c(1,2,3,4,5,6,7,8,9,10),"Feature_Name"] # keeps only the genes on one of the ten chromosomes

save_random=T
#Duplicate taxa avearged or not averaged?
dup_avg=F
#All >0 median exp genes and top 5000 genes, and 5001-10000, and busco, and random 5000 genes
for (j in unique(expression_df$TissueWODate)){
  cat(j, " ", as.character(Sys.time()), "\n")
  # pull genes from gff file and keep only the genes which are on chr1 - chr10
  # gff_file=read.delim("/home/kak268/work_data/maize_genome/Zea_mays_gene_only.AGPv3.29.gff3", header=F, sep="\t", quote="") # add skip=6 if it's a full gff file
  # colnames(gff_file)=c("Chr", "Annot_source", "Feature_type", "Start", "End", "blankcol1", "Strand", "blankcol2", "Feature_Name")
  # transcript_file=gff_file[gff_file$Feature_type=='gene',] # keep only the entries in the GFF file which are transcripts/genes
  # transcript_file$Feature_Name=sub(pattern= "ID=gene:", replacement='', x=transcript_file$Feature_Name)
  # transcript_file$Feature_Name=sub(pattern= ";.*", replacement='', x=transcript_file$Feature_Name)
  # transcript_file_only_chr1_to_chr10_gene_list=transcript_file[transcript_file$Chr %in% c(1,2,3,4,5,6,7,8,9,10),"Feature_Name"] # keeps only the genes on one of the ten chromosomes
  
  
  SingleTissue_df=expression_df[expression_df$TissueWODate==j,]
  #write.table(SinglTissue_df,paste(basedir,"L3tip_vs_rare_alleles/df_STAR_HTSeq_counts_B73_match_based_on_genet_dist_DESeq2_normed_rounded_origNames_and_Altnames_L3Tip.txt", sep=""), quote=F, row.names=F)
  SingleTissue_df[1:3,1:15]
  SingleTissue_df_w_taxa_only=SingleTissue_df[-c(1:5, 7:10)]
  
  dup_avg_string_for_filenames="DupTaxaNotAveraged"
  if (dup_avg==T){ #this allows you to turn avg over dups on and off
    dup_avg_string_for_filenames="AvgOverDups"
    #this averages over all duplicate taxa, IN PARALLEL!!!
    SingleTissue_df_avg_over_taxa=ddply(SingleTissue_df_w_taxa_only, "HMP32Name", numcolwise(mean), .parallel = T)
    SingleTissue_df_w_taxa_only=SingleTissue_df_avg_over_taxa
  }
  
  #Strips off "HMP32Name" col and keeps only cols which are for genes on chr1-10
  Non_redundant_taxa_names=SingleTissue_df_w_taxa_only$HMP32Name
  SingleTissue_df_no_names=SingleTissue_df_w_taxa_only[,colnames(SingleTissue_df_w_taxa_only) %in% transcript_file_only_chr1_to_chr10_gene_list]
  #Remove the 
  #expression_df[expression_df=="NaN"]=NA
  
  gene_name_with_median_expr=as.data.frame(cbind(Feature_Name=colnames(SingleTissue_df_no_names), expr_median=colMedians(SingleTissue_df_no_names)), stringsAsFactors=F)
  gene_name_with_median_expr_no_0=gene_name_with_median_expr[gene_name_with_median_expr$expr_median!=0,]
  sorted_gene_name_with_median_expr=gene_name_with_median_expr_no_0[order(-(as.numeric(gene_name_with_median_expr_no_0[,2]))),]
  #hist(log(as.numeric(sorted_gene_name_with_median_expr$expr_median)))
  
  #retrieve all non0 median expressed
  all_non0_genes=sorted_gene_name_with_median_expr
  write.table(all_non0_genes$Feature_Name, paste(basedir, "all_non0_median_exp_genes_by_tissue_",dup_avg_string_for_filenames,"/", j, "_all_non0_median_STAR_HTSEQ_DESEQ_expressed_genes.txt", sep=""), quote=F, row.names=T, col.names=T)
  write.table(all_non0_genes, paste(basedir, "all_non0_median_exp_genes_by_tissue_",dup_avg_string_for_filenames,"/", j, "_all_non0_expressed_genes_w_STAR_HTSEQ_DESEQ_median_counts.txt", sep=""), quote=F, row.names=F, col.names=T)
  #Write the entire expression matrix for these to file
  SingleTissue_df_avg_over_taxa_w_nonRedTaxaNames=cbind(Non_redundant_taxa_names,SingleTissue_df_w_taxa_only)
  all_non0_gene_w_exr_values=SingleTissue_df_avg_over_taxa_w_nonRedTaxaNames[, colnames(SingleTissue_df_avg_over_taxa_w_nonRedTaxaNames) %in% c("HMP32Name", all_non0_genes$Feature_Name)]
  write.table(all_non0_gene_w_exr_values, paste(basedir, "all_non0_median_exp_genes_by_tissue_",dup_avg_string_for_filenames,"/", j, "_all_non0_expressed_genes_w_STAR_HTSEQ_DESEQ_actual_counts.txt", sep=""), sep="\t", quote=F, row.names=F, col.names=T)
  
  
  
  #retrieve top 5000 median expressed genes and save values
  top_5000_genes=sorted_gene_name_with_median_expr[1:5000,]
  write.table(top_5000_genes$Feature_Name, paste(basedir, "Top_5000_by_tissue_",dup_avg_string_for_filenames,"/", j, "_Top_5000_median_STAR_HTSEQ_DESEQ_expressed_genes.txt", sep=""), quote=F, row.names=T, col.names=T)
  write.table(top_5000_genes, paste(basedir, "Top_5000_by_tissue_",dup_avg_string_for_filenames,"/", j, "_Top_5000_expressed_genes_w_STAR_HTSEQ_DESEQ_median_counts.txt", sep=""), quote=F, row.names=F, col.names=T)
  #Write the entire expression matrix for these top 5000 to file
  SingleTissue_df_avg_over_taxa_w_nonRedTaxaNames=cbind(Non_redundant_taxa_names,SingleTissue_df_w_taxa_only)
  top_5000_gene_w_exr_values=SingleTissue_df_avg_over_taxa_w_nonRedTaxaNames[, colnames(SingleTissue_df_avg_over_taxa_w_nonRedTaxaNames) %in% c("HMP32Name", top_5000_genes$Feature_Name)]
  write.table(top_5000_gene_w_exr_values, paste(basedir, "Top_5000_by_tissue_",dup_avg_string_for_filenames,"/", j, "_Top_5000_expressed_genes_w_STAR_HTSEQ_DESEQ_actual_counts.txt", sep=""), sep="\t", quote=F, row.names=F, col.names=T)
  # 
  #retrieve top 5001_10000 median expressed genes and save values
  top_5001_to_10000_genes=sorted_gene_name_with_median_expr[5001:10000,]
  write.table(top_5001_to_10000_genes$Feature_Name, paste(basedir, "Top_5001_to_10000_by_tissue_",dup_avg_string_for_filenames,"/", j, "_Top_5001_to_10000_median_STAR_HTSEQ_DESEQ_expressed_genes.txt", sep=""), quote=F, row.names=T, col.names=T)
  write.table(top_5001_to_10000_genes, paste(basedir, "Top_5001_to_10000_by_tissue_",dup_avg_string_for_filenames,"/", j, "_Top_5001_to_10000_expressed_genes_w_STAR_HTSEQ_DESEQ_median_counts.txt", sep=""), quote=F, row.names=F, col.names=T)
  #Write the entire expression matrix for these top 5001_10000 to file
  SingleTissue_df_avg_over_taxa_w_nonRedTaxaNames=cbind(Non_redundant_taxa_names,SingleTissue_df_w_taxa_only)
  top_5001_to_10000_gene_w_exr_values=SingleTissue_df_avg_over_taxa_w_nonRedTaxaNames[, colnames(SingleTissue_df_avg_over_taxa_w_nonRedTaxaNames) %in% c("HMP32Name", top_5001_to_10000_genes$Feature_Name)]
  write.table(top_5001_to_10000_gene_w_exr_values, paste(basedir, "Top_5001_to_10000_by_tissue_",dup_avg_string_for_filenames,"/", j, "_Top_5001_to_10000_expressed_genes_w_STAR_HTSEQ_DESEQ_actual_counts.txt", sep=""), sep="\t", quote=F, row.names=F, col.names=T)
  
  
  #retrieve BUSCO geens with median expression over 0
  busco_genes=read.table(file = "/media/kak268/B_2TB_Internal/BUSCO/BUSCO_Zea_Mays.AGPv3.29.pep.all/run_out_BUSCO_Zea_Mays.AGPv3.29.pep.all/Unique_genes_full_table_out_BUSCO_Zea_Mays.AGPv3.29.pep.all", header=F)
  busco_genes_overlap=sorted_gene_name_with_median_expr[rownames(sorted_gene_name_with_median_expr) %in% busco_genes[,1],]
  write.table(busco_genes_overlap$Feature_Name, paste(basedir, "busco_by_tissue_",dup_avg_string_for_filenames,"/", j, "_busco_non0_median_STAR_HTSEQ_DESEQ_expressed_genes.txt", sep=""), quote=F, row.names=T, col.names=T)
  write.table(busco_genes_overlap, paste(basedir, "busco_by_tissue_",dup_avg_string_for_filenames,"/", j, "_busco_non0_expressed_genes_w_STAR_HTSEQ_DESEQ_median_counts.txt", sep=""), quote=F, row.names=F, col.names=T)
  #Write the entire expression matrix for these top 5001_10000 to file
  SingleTissue_df_avg_over_taxa_w_nonRedTaxaNames=cbind(Non_redundant_taxa_names,SingleTissue_df_w_taxa_only)
  busco_gene_w_exr_values=SingleTissue_df_avg_over_taxa_w_nonRedTaxaNames[, colnames(SingleTissue_df_avg_over_taxa_w_nonRedTaxaNames) %in% c("HMP32Name", busco_genes_overlap$Feature_Name)]
  write.table(busco_gene_w_exr_values, paste(basedir, "busco_by_tissue_",dup_avg_string_for_filenames,"/", j, "_busco_non0_expressed_genes_w_STAR_HTSEQ_DESEQ_actual_counts.txt", sep=""), sep="\t", quote=F, row.names=F, col.names=T)
  
  
  
  #Retrieve random 5000 genes and save values
  if (save_random==T){
    set.seed(123)
    rand_5000_non0_genes=sorted_gene_name_with_median_expr[sample(nrow(sorted_gene_name_with_median_expr), 5000),]
    write.table(rand_5000_non0_genes$Feature_Name, paste(basedir, "Rand_5000_by_tissue_",dup_avg_string_for_filenames,"/", j, "_Rand_5000_median_STAR_HTSEQ_DESEQ_expressed_genes.txt", sep=""), quote=F, row.names=T, col.names=T)
    write.table(rand_5000_non0_genes, paste(basedir, "Rand_5000_by_tissue_",dup_avg_string_for_filenames,"/", j, "_Rand_5000_expressed_genes_w_STAR_HTSEQ_DESEQ_median_counts.txt", sep=""), quote=F, row.names=F, col.names=T)
    #Write the entire expression matrix for these rand 5000 to file
    SingleTissue_df_avg_over_taxa_w_nonRedTaxaNames=cbind(Non_redundant_taxa_names,SingleTissue_df_w_taxa_only)
    rand_5000_gene_w_exr_values=SingleTissue_df_avg_over_taxa_w_nonRedTaxaNames[, colnames(SingleTissue_df_avg_over_taxa_w_nonRedTaxaNames) %in% c("HMP32Name", rand_5000_non0_genes$Feature_Name)]
    write.table(rand_5000_gene_w_exr_values, paste(basedir, "Rand_5000_by_tissue_",dup_avg_string_for_filenames,"/", j, "_Rand_5000_expressed_genes_w_STAR_HTSEQ_DESEQ_actual_counts.txt", sep=""), sep="\t", quote=F, row.names=F, col.names=T)
  }
  
}

