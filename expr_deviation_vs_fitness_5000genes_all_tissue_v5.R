#Karl Kremling
#Correlate mean(abs( val of expr rank - mean expr rank)) vs fitness (height or yield)

library(reshape)
library(miscTools)
library(plyr)
library("foreach")
library("doMC")
#library("parallel")
registerDoMC(cores=12)
library(reshape2)
library(readr)
#install.packages("heR.misc") # for geometric mean calculations
#library(heR.misc)

#Calculate the expression RESIDUALS after factoring out Hidden Factors calculated from expression

tissue_list=c("Kern", "GRoot", "L3Base", "L3Tip", "LMAD", "LMAN", "GShoot")
#tissue_list=c("Kern", "GRoot", "L3Base")



expression_df=read.table("/media/kak268/A_2TB_Internal/RNA/Expressions_quants/HTSEQ_counts_STAR_against_Zea_mays_B73_AGPv3.29_w_gtf_annotation/count_mat/count_mats_matched_to_AltNames/df_STAR_HTSeq_counts_B73_match_based_on_genet_dist_DESeq2_normed_rounded_origNames_and_Altnames_noL3Mid.txt", header = T)


#Test head2=head(PEER_df_no_MDS[,1:6])
#Test head1=head(single_tissue_expression[,1:6])

###
###
#Calculate residual after removing PEERs
#This uses PEER residuals and factors them out. 
top_5000_only=T
for(h in tissue_list){
  rm()
  rm(res)
  PEER_df=read.table(paste("/media/kak268/B_2TB_Internal/Genotypes/HM32/LDKNNi/PCs/",h,"/MDSPCs_from_HMP321_w_25_PEER_factors_",h,"_complete_taxa_and_tiss_names.txt",sep = ""),header = T)
  PEER_df_no_MDS=PEER_df[,c(1,16:dim(PEER_df)[2])]
  #single_tissue_expression=expression_df[grep(pattern = as.character(h), expression_df$X.Trait.),]
  single_tissue_expression=expression_df[grep(pattern = as.character(h), expression_df$TissueWODate),]
  if(top_5000_only==T){
    top_5000_list=read.table(paste("/media/kak268/A_2TB_Internal/RNA/Expressions_quants/HTSEQ_counts_STAR_against_Zea_mays_B73_AGPv3.29_w_gtf_annotation/count_mat/count_mats_matched_to_AltNames/Top_5000_by_tissue_DupTaxaNotAveraged/", h, "_Top_5000_expressed_genes_w_STAR_HTSEQ_DESEQ_actual_counts.txt", sep = ""),nrows = 1, header = F)
  
    col_idx_top_5000_list=which(colnames(single_tissue_expression) %in% t(top_5000_list[1,-1])) # you MUST transpose the df first and skip repeated "HMP32Name"
    
    #Keep only the top 5000 columns and stick the unique identifier back on
    single_tissue_expression_top_5000=cbind(OrigColNameFromRNAExpressionValue_RNA_TaxaName=single_tissue_expression$OrigColNameFromRNAExpressionValue_RNA_TaxaName,single_tissue_expression[,col_idx_top_5000_list])
    
    #Merge expression and PEERs then separate the dfs to make sure they are in the same sorted order    
    merged_df=merge(single_tissue_expression_top_5000, PEER_df_no_MDS, by= "OrigColNameFromRNAExpressionValue_RNA_TaxaName")
  }
  if(top_5000_only==F){
    
    non0_median_list=read.table(paste("/media/kak268/A_2TB_Internal/RNA/Expressions_quants/HTSEQ_counts_STAR_against_Zea_mays_B73_AGPv3.29_w_gtf_annotation/count_mat/count_mats_matched_to_AltNames/all_non0_median_exp_genes_by_tissue_DupTaxaNotAveraged/", h, "_all_non0_expressed_genes_w_STAR_HTSEQ_DESEQ_actual_counts.txt", sep = ""),nrows = 1, header = F)
    
    col_idx_non0_median_list=which(colnames(single_tissue_expression) %in% t(non0_median_list[1,-1])) # you MUST transpose the df first and skip repeated "HMP32Name"
    #Keep only the top 5000 columns and stick the unique identifier back on
    single_tissue_expression_all_non0=cbind(OrigColNameFromRNAExpressionValue_RNA_TaxaName=single_tissue_expression$OrigColNameFromRNAExpressionValue_RNA_TaxaName, single_tissue_expression[,col_idx_non0_median_list])
    
    #Merge expression and PEERs then separate the dfs to make sure they are in the same sorted order    
    merged_df=merge(single_tissue_expression_all_non0, PEER_df_no_MDS, by= "OrigColNameFromRNAExpressionValue_RNA_TaxaName")
    
  }
  rownames(merged_df)=merged_df$OrigColNameFromRNAExpressionValue_RNA_TaxaName
  
  merged_df_expr=merged_df[,2:(dim(merged_df)[2]-(dim(PEER_df_no_MDS)[2]-1))]
  merged_df_PEERs=merged_df[,(dim(merged_df)[2]-(dim(PEER_df_no_MDS)[2]-2)):(dim(merged_df)[2])]
  
  #calculate residuals
  res=residuals(lm(as.matrix(merged_df_expr)~as.matrix(merged_df_PEERs)))
  # convert rownames into first column and cbind HMP32 name again
  res = cbind.data.frame(OrigColNameFromRNAExpressionValue_RNA_TaxaName=rownames(res), res)
  res=merge(single_tissue_expression[,c("OrigColNameFromRNAExpressionValue_RNA_TaxaName", "HMP32Name")],res)
  assign(paste("residuals", h, sep=""), res)
  #Re-add the HMP32 names to the residual df
  if(top_5000_only==T){
    write.table(x = res, quote = F, row.names = F, file =  paste("/media/kak268/A_2TB_Internal/RNA/Expressions_quants/HTSEQ_counts_STAR_against_Zea_mays_B73_AGPv3.29_w_gtf_annotation/count_mat/count_mats_matched_to_AltNames/Top_5000_by_tissue_DupTaxaNotAveraged_residuals_after_removing_PEERs/",h,"_Top_5000_expressed_genes_w_STAR_HTSEQ_DESEQ_residuals_after_25PEERs_factored_out.txt", sep = ""))
  }
  if(top_5000_only==F){
    write.table(x = res, quote = F, row.names = F, file =  paste("/media/kak268/A_2TB_Internal/RNA/Expressions_quants/HTSEQ_counts_STAR_against_Zea_mays_B73_AGPv3.29_w_gtf_annotation/count_mat/count_mats_matched_to_AltNames/all_non0_median_genes_residuals_after_removing_PEERs_DupTaxaNotAveraged/",h,"_all_non0_median_w_STAR_HTSEQ_DESEQ_residuals_after_25PEERs_factored_out.txt", sep = ""))
  }
}


#Top 5000 genes from each tissue, dups not averaged
basedir="/media/kak268/A_2TB_Internal/RNA/Expressions_quants/HTSEQ_counts_STAR_against_Zea_mays_B73_AGPv3.29_w_gtf_annotation/count_mat/count_mats_matched_to_AltNames/Top_5000_by_tissue_DupTaxaNotAveraged/"
filesuffix="in 5000 most exp genes 0notSetToNA_DupTaxaNotAveraged"

#this is what is used in the paper
#Top 5000 genes from each tissue avg over dups
basedir="/media/kak268/A_2TB_Internal/RNA/Expressions_quants/HTSEQ_counts_STAR_against_Zea_mays_B73_AGPv3.29_w_gtf_annotation/count_mat/count_mats_matched_to_AltNames/Top_5000_by_tissue_AvgOverDups/"
filelist=list.files(basedir, pattern="actual_counts.txt")
plot_out_stemdir="/home/kak268/work_data/maize_phenotypes/Yield_BLUPs_on_Testers/plots/expression_dev_v_top_5000_w_sweet_and_pop_exclusions"
exclude=T
#exclusion_list=c("282set_Ia5125", "282set_Il101", "282set_Il14H", "282set_P39Goodman-Buckler", "282set_i1677a", "282set_IA2132Goodman-Buckler")
exclusion_list=c("282set_Ia5125", "282set_Il101", "282set_Il14H", "282set_P39Goodman-Buckler", "282set_i1677a", "282set_IA2132Goodman-Buckler", "282set_IDS91", "282set_4722", "282set_HP301", "282set_I29", "282set_IDS28", "282set_IDS69", "282set_IDS91", "282set_SA24", "282set_Sg1533", "282set_Sg18")
filesuffix="in 5k most exp genes 0notSetToNA_AvgOverDups no sweet no pop"
only_over_or_under_exp=F


#Top 5000 genes from each tissue avg over dups, SEPARATE over and under expression
basedir="/media/kak268/A_2TB_Internal/RNA/Expressions_quants/HTSEQ_counts_STAR_against_Zea_mays_B73_AGPv3.29_w_gtf_annotation/count_mat/count_mats_matched_to_AltNames/Top_5000_by_tissue_AvgOverDups/"
filelist=list.files(basedir, pattern="actual_counts.txt")
plot_out_stemdir="/home/kak268/work_data/maize_phenotypes/Yield_BLUPs_on_Testers/plots/expression_dev_v_top_5000_w_sweet_and_pop_exclusions_only_over_or_under_exp"
exclude=T
#exclusion_list=c("282set_Ia5125", "282set_Il101", "282set_Il14H", "282set_P39Goodman-Buckler", "282set_i1677a", "282set_IA2132Goodman-Buckler")
exclusion_list=c("282set_Ia5125", "282set_Il101", "282set_Il14H", "282set_P39Goodman-Buckler", "282set_i1677a", "282set_IA2132Goodman-Buckler", "282set_IDS91", "282set_4722", "282set_HP301", "282set_I29", "282set_IDS28", "282set_IDS69", "282set_IDS91", "282set_SA24", "282set_Sg1533", "282set_Sg18")
only_over_or_under_exp=F
pos_or_neg="pos"
filesuffix=paste("in 5k most exp genes 0notSetToNA_AvgOverDups no sweet no pop only ", pos_or_neg, sep="")



#All non 0 genes from each tissue
basedir="/media/kak268/A_2TB_Internal/RNA/Expressions_quants/HTSEQ_counts_STAR_against_Zea_mays_B73_AGPv3.29_w_gtf_annotation/count_mat/count_mats_matched_to_AltNames/all_non0_median_exp_genes_by_tissue_DupTaxaNotAveraged/"
filelist=list.files(basedir, pattern="actual_counts.txt")
plot_out_stemdir="/home/kak268/work_data/maize_phenotypes/Yield_BLUPs_on_Testers/plots/expression_dev_all_non_0_0notSetToNA"
filesuffix="in all median > 0 exp genes 0notSetToNA"

#Top 5000 genes from each tissue after regressing out PEER factors, dups not YET averaged (dups averaged in this script)
basedir="/media/kak268/A_2TB_Internal/RNA/Expressions_quants/HTSEQ_counts_STAR_against_Zea_mays_B73_AGPv3.29_w_gtf_annotation/count_mat/count_mats_matched_to_AltNames/Top_5000_by_tissue_DupTaxaNotAveraged_residuals_after_removing_PEERs/"
filelist=list.files(basedir, pattern="residuals_after_25PEERs_factored_out.txt")
filesuffix="in 5000 most exp genes 0notSetToNA"
plot_out_stemdir="/home/kak268/work_data/maize_phenotypes/Yield_BLUPs_on_Testers/plots/expression_dev_v_top_5000_0notSetToNA_using_residuals"


#All non0 median genes from each tissue after regressing out PEER factors, dups not YET averaged (dups averaged in this script)
basedir="/media/kak268/A_2TB_Internal/RNA/Expressions_quants/HTSEQ_counts_STAR_against_Zea_mays_B73_AGPv3.29_w_gtf_annotation/count_mat/count_mats_matched_to_AltNames/all_non0_median_genes_residuals_after_removing_PEERs_DupTaxaNotAveraged/"
filelist=list.files(basedir, pattern="residuals_after_25PEERs_factored_out.txt")
filesuffix="in all median > 0 exp genes 0notSetToNA"
plot_out_stemdir="/home/kak268/work_data/maize_phenotypes/Yield_BLUPs_on_Testers/plots/expression_dev_v_all_non0_median_0notSetToNA_using_residuals"


#plot_out_dir="/home/kak268/work_data/maize_phenotypes/Yield_BLUPs_on_Testers/plots/expression_dev_v_top_5000/"
#plot_out_dir="/home/kak268/work_data/maize_phenotypes/Yield_BLUPs_on_Testers/plots/expression_dev_all_non_0/"



use_residuals=F
avg_over_dups=T
if (avg_over_dups==T){ #only needed if data haven't already been averaged over duplicates
  plot_out_dir=paste(plot_out_stemdir, "_AvgOverDups/", sep="")
  cat("avg over dups status is ", avg_over_dups)
}else{
  plot_out_dir=paste(plot_out_stemdir, "_DupTaxaNotAveraged/", sep="")
}
for (i in filelist){
  cat(i, "\n")
  print(Sys.time())
  if(exists("gene_w_exr_values")){
    rm(gene_w_exr_values)
  }
  if (i!="L3Mid_top_5000_expressed_genes_w_STAR_HTSEQ_DESEQ_actual_counts.txt"){
    tissue_name=substr(i, 1, which(strsplit(i, "")[[1]]=="_")[1]-1)
    if(use_residuals==F){
      gene_w_exr_values=read.table(paste(basedir, i, sep=""), header=T)
      cat("gene w exr values before avg before exclusion list applied \n")
      cat(dim(gene_w_exr_values))
      if (exclude==T){
        gene_w_exr_values=gene_w_exr_values[! gene_w_exr_values$HMP32Name %in% exclusion_list,]
      }
      cat("gene w exr values before avg after exclusion list applied \n")
      cat(dim(gene_w_exr_values))
      if (avg_over_dups==T){
        SingleTissue_df_avg_over_taxa=ddply(gene_w_exr_values, "HMP32Name", numcolwise(mean), .parallel = T)
        gene_w_exr_values=SingleTissue_df_avg_over_taxa
        cat("gene w exr values after avg")
        cat(dim(gene_w_exr_values))
      }
    }
    if(use_residuals==T){
      gene_w_exr_values=read.table(file = paste(basedir, i, sep=""), header=T)
      gene_w_exr_values=gene_w_exr_values[,-1]
      if (avg_over_dups==T){
        SingleTissue_df_avg_over_taxa=ddply(gene_w_exr_values, "HMP32Name", numcolwise(mean), .parallel = T)
        gene_w_exr_values=SingleTissue_df_avg_over_taxa
      }
    }
    
    #Pre calc the deviations
    means=colMeans(gene_w_exr_values[,2:ncol(gene_w_exr_values)])
    var_mat=gene_w_exr_values
    var_mat[,-1]=NA # make an empty mat of same dimensions, leave names in first column
    for (j in 1:nrow(gene_w_exr_values)){
      if(only_over_or_under_exp==F){
      var_mat[j,2:ncol(var_mat)]=(gene_w_exr_values[j,2:ncol(var_mat)]-means)^2
      }else if(only_over_or_under_exp==T){
        if(pos_or_neg=="pos"){
        var_mat[j,2:ncol(var_mat)]=gene_w_exr_values[j,2:ncol(var_mat)]-means
        #var_mat[j,2:ncol(var_mat)]<0=0
        var_mat[j,c(FALSE,var_mat[j,2:ncol(var_mat)]<0)]=0
        var_mat[j,2:ncol(var_mat)]=(var_mat[j,2:ncol(var_mat)])^2
        }else if(pos_or_neg=="neg"){
          var_mat[j,2:ncol(var_mat)]=gene_w_exr_values[j,2:ncol(var_mat)]-means
          #var_mat[j,2:ncol(var_mat)]>0=0
          var_mat[j,c(FALSE,var_mat[j,2:ncol(var_mat)]>0)]=0
          var_mat[j,2:ncol(var_mat)]=(var_mat[j,2:ncol(var_mat)])^2
        }
      }
      cat(j)
    }
    
    plot_exp_dev_v_phenos_mean(gene_w_exr_values=gene_w_exr_values, tissue=tissue_name, xlab=paste("log10 Mean sq. deviation of individual ", tissue_name, " expression from population mean ", filesuffix, sep=""), plot_out_dir)
    #plot_exp_dev_v_phenos_median(gene_w_exr_values=gene_w_exr_values, tissue=tissue_name, xlab=paste("log10 Median sq. deviation of individual ", tissue_name, " expression from population mean ", filesuffix, sep=""), plot_out_dir) #" expression from population mean ", filesuffix, sep=""))
    if(filesuffix!="in all median > 0 exp genes 0notSetToNA"){ # you cannot do geom means if too many zeroes in the exp values without setting them to NA
      #plot_exp_dev_v_phenos_geo_mean(gene_w_exr_values=gene_w_exr_values, tissue=tissue_name, xlab=paste("log10 Geom mean sq. deviation of individual ", tissue_name, " expression from population mean ", filesuffix, sep=""), plot_out_dir) #" expression from population mean ", filesuffix, sep=""))
    }
    
    cat("\n", "Starting rank based correl", "\n")
    
    #RANKBASED Pre calc the deviations in rank
    gene_w_exr_values_ranks=gene_w_exr_values
    gene_w_exr_values_ranks[,-1]=NA
    gene_w_exr_values_ranks[,2:ncol(gene_w_exr_values_ranks)]=apply(gene_w_exr_values[,2:ncol(gene_w_exr_values)], 2,function(x) rank(x))
    means=colMeans(gene_w_exr_values_ranks[,2:ncol(gene_w_exr_values_ranks)])
    var_mat_rank=gene_w_exr_values
    var_mat_rank[,-1]=NA
    for (i in 1:nrow(gene_w_exr_values)){
      var_mat_rank[i,2:ncol(var_mat_rank)]=(gene_w_exr_values[i,2:ncol(var_mat_rank)]-means)^2
      cat(i)
    }
    #plot_exp_dev_v_phenos_rank(gene_w_exr_values=gene_w_exr_values, tissue=tissue_name, xlab=paste("log10 Mean sq. deviation of individual ", tissue_name, " expression from population mean ", filesuffix, sep=""), plot_out_dir) #" rank expression from population mean ", filesuffix, sep=""))
  }
}




print(Sys.time())
#Use all tissues together for plotting
plot_exp_dev_v_phenos_mean(gene_w_exr_values=dev_mat_master_mean, tissue="all_tissues", xlab=paste("log10 Mean sq. deviation of individual ", "all_tissues", " expression from population mean ", filesuffix, sep=""), plot_out_dir)
plot_exp_dev_v_phenos_median(gene_w_exr_values=dev_mat_master_median, tissue="all_tissues", xlab=paste("log10 Median sq. deviation of individual ", "all_tissues", " expression from population mean ", filesuffix, sep=""), plot_out_dir)
if(filesuffix!="in all median > 0 exp genes 0notSetToNA"){ # you cannot do geom means if too many zeroes in the exp values without setting them to NA
  plot_exp_dev_v_phenos_geo_mean(gene_w_exr_values=dev_mat_master_geo_mean, tissue="all_tissues", xlab=paste("log10 Geom mean sq. deviation of individual ", "all_tissues", " expression from population mean ", filesuffix, sep=""), plot_out_dir)
}
plot_exp_dev_v_phenos_rank(gene_w_exr_values=dev_mat_rank_master_mean, tissue="all_tissues", xlab=paste("log10 Mean sq. deviation of individual rank ", "all_tissues", " expression from population mean ", filesuffix, sep=""), plot_out_dir)
print(Sys.time())


##These are the plotting functions
#Expr value based
#with mean
dev_mat_master_mean=NULL
plot_exp_dev_v_phenos_mean=function(gene_w_exr_values, tissue, xlab, plot_out_dir){
  dev_mat_master_mean=NULL
  if(tissue!="all_tissues"){
    #Mean dev over all genes in each taxon
    #MEAN
    #mean_dev_by_taxon=cbind.data.frame(HMP32Name=as.character(var_mat$HMP32Name), Mean_Dev=as.numeric(rowMeans(var_mat[,2:ncol(var_mat)])))
    #MEDIAN
    #mean_dev_by_taxon=cbind.data.frame(HMP32Name=as.character(var_mat$HMP32Name), Mean_Dev=as.numeric(rowMedians(var_mat[,2:ncol(var_mat)])))
    #GEOMETRIC MEAN
    #mean_dev_by_taxon=cbind.data.frame(HMP32Name=as.character(var_mat$HMP32Name), Mean_Dev=as.numeric(apply(var_mat[,2:ncol(var_mat)], 1, function(x) exp(mean(log(x), na.rm=T)))))
    
    
    #Remove genes w 0 exp
    #NOTE THE ZERO EXP GENE REMOVAL IS NOW OMITTED BECAUSE WE'RE JUST USING TOP 5000 genes, BUT I KEPT THE VARIABLE NAME
    #To ensure that distantly related individuals which simply lack some B73 genes don't have inflated deviations, remove the variances from those genes which have zero expression
    var_mat_zero_exp_genes_set_to_zero_var=var_mat
    #var_mat_zero_exp_genes_set_to_zero_var[gene_w_exr_values==0]=NA
  
    #MEAN
    mean_dev_by_taxon_exc_zero_exp_genes=cbind.data.frame(HMP32Name=as.character(var_mat_zero_exp_genes_set_to_zero_var$HMP32Name), Mean_Dev=as.numeric(rowMeans(var_mat_zero_exp_genes_set_to_zero_var[,2:ncol(var_mat_zero_exp_genes_set_to_zero_var)], na.rm=T)))
    #MEDIAN
    #mean_dev_by_taxon_exc_zero_exp_genes=cbind.data.frame(HMP32Name=as.character(var_mat_zero_exp_genes_set_to_zero_var$HMP32Name), Mean_Dev=as.numeric(rowMedians(var_mat_zero_exp_genes_set_to_zero_var[,2:ncol(var_mat_zero_exp_genes_set_to_zero_var)], na.rm=T)))
    #GEOMETRIC MEAN
    #mean_dev_by_taxon_exc_zero_exp_genes=cbind.data.frame(HMP32Name=as.character(var_mat_zero_exp_genes_set_to_zero_var$HMP32Name), Mean_Dev=as.numeric(apply(var_mat_zero_exp_genes_set_to_zero_var[,2:ncol(var_mat_zero_exp_genes_set_to_zero_var)], 1, function(x) exp(mean(log(x), na.rm=T)))))
    
    #Seed Weight BLUPs
    #Load in phenotypes 
    #phenos=read.table("/home/kak268/work_data/maize_phenotypes/Ames_BLUPs/Merged_Ames_traits", header=T)
    #phenos=read.table(file="/home/kak268/work_data/maize_phenotypes/282_panzea_phenotypes/traitMatrix_maize282NAM_v15-130212_Earweight_minus_Cobweight_melt_Means_and_BLUPS.txt", sep="\t", header=T)
    
    
    #Load in name converter so names can be converted to Ames names
    name_decoder=read.table(header=T,"/media/kak268/A_2TB_Internal/RNA/Expressions_quants/renaming_scratch_folder/orig_1960names_and_alternate_names.txt")  
    name_decoder_tissue=name_decoder[name_decoder$TissueWODate==tissue,]
    ###Make the decoder unique
    unique_name_decoder_tissue=name_decoder_tissue[!duplicated(name_decoder_tissue$HMP32Name),]
    
    
    #Merge dev by taxon with name decoder and phenos
    #inc zero exp gene dev
    #mean_dev_by_taxon_w_names=merge(mean_dev_by_taxon, unique_name_decoder_tissue, by="HMP32Name")
    #exc zero exp gene dev
    mean_dev_by_taxon_w_names=merge(mean_dev_by_taxon_exc_zero_exp_genes, unique_name_decoder_tissue, by="HMP32Name")
    
    mean_dev_by_taxon_w_names_unique_combs=mean_dev_by_taxon_w_names[!duplicated(mean_dev_by_taxon_w_names[c("HMP32Name", "AmesName", "Mean_Dev", "RNA_TaxaName", "RawPhenotypeNames")]),]
    dev_mat_master_mean=rbind(dev_mat_master_mean, mean_dev_by_taxon_w_names_unique_combs)
    assign('dev_mat_master_mean', dev_mat_master_mean, envir=.GlobalEnv)
  } else{ # if tissue=="all_tissues"
    mean_dev_by_taxon_w_names_unique_combs=gene_w_exr_values    # in this case this holds the precomputed deviations from each tissue
  }
  
  #NO longer using the means fo the BLUPS in each enviro
  ##Seed Weight Means
  ##Load in phenotypes 
  # phenos=read.table(file="/home/kak268/work_data/maize_phenotypes/282_panzea_phenotypes/traitMatrix_maize282NAM_v15-130212_Earweight_minus_Cobweight_melt_Means_and_BLUPS.txt", sep="\t", header=T)
  # 
  # #orig#mean_dev_by_taxon_w_names_and_phenos_unique_combs=merge(x=mean_dev_by_taxon_w_names_unique_combs, y=phenos, by.x="AmesName", by.y="TaxaName") 
  # mean_dev_by_taxon_w_names_and_phenos_unique_combs=merge(x=mean_dev_by_taxon_w_names_unique_combs, y=phenos, by.x="RawPhenotypeNames", by.y="TaxaName") 
  # 
  #   #orig plot(mean_dev_by_taxon_w_names_and_phenos_unique_combs$Mean_Dev, mean_dev_by_taxon_w_names_and_phenos_unique_combs$PHT_BLUP, cex=1.5, col="deeppink4", pch=16, xlab=xlab, ylab=ylab)
  # pdf(paste(plot_out_dir,xlab,"vs_mean_seed_wt_BLUP.pdf", sep=""),width=11,height=11,paper='special') 
  # plot(log10(mean_dev_by_taxon_w_names_and_phenos_unique_combs$Mean_Dev), mean_dev_by_taxon_w_names_and_phenos_unique_combs$SeedWeightBLUPs, cex=1.5, col="darkorange3", pch=16, xlab=xlab, ylab="Seed Weight BLUP means")
  # #orig line=lm(mean_dev_by_taxon_w_names_and_phenos_unique_combs$PHT_BLUP~mean_dev_by_taxon_w_names_and_phenos_unique_combs$Mean_Dev)
  # line=lm(mean_dev_by_taxon_w_names_and_phenos_unique_combs$SeedWeightBLUPs~log10(mean_dev_by_taxon_w_names_and_phenos_unique_combs$Mean_Dev))
  # abline(line)
  # linesummary=summary(line)
  # #legend("topright", bty="n", legend=as.expression(c(bquote(R^2 == .(round(linesummary$r.squared, digits=4))), paste("P = ", bquote(.(round(pf(linesummary$fstatistic[1], linesummary$fstatistic[2], linesummary$fstatistic[3], lower.tail = FALSE), digits=10))), sep=""))))# p == .(pval1) .(pval2)))
  # legend("topright", bty="n", legend=as.expression(c(bquote(r == .(round(sqrt(linesummary$r.squared), digits=4))), paste("P = ", bquote(.(round(pf(linesummary$fstatistic[1], linesummary$fstatistic[2], linesummary$fstatistic[3], lower.tail = FALSE), digits=10))), sep=""))))# p == .(pval1) .(pval2)))
  # 
  # intercept=round(linesummary$coefficients[1,1], 2)
  # slope=round(linesummary$coefficients[2,1], 2)
  # sign=NULL
  # if (slope>=0){sign="+"}
  # legend("top", bty="n", legend=paste("y = ", intercept, sign, slope,"x", sep=""))
  # dev.off()
  
  #HOLLAND 06 07 Seed Weight BLUPs
  rm(phenos)
  #Load in phenotypes 
  #phenos=read.table("/home/kak268/work_data/maize_phenotypes/Ames_BLUPs/Merged_Ames_traits", header=T)
  phenos=read.table(file="/home/kak268/work_data/maize_phenotypes/Holland_NAM_and_282_BLUPs/282 BLUPs 0607 all traits.txt", sep="\t", header=T)
  
  #orig#mean_dev_by_taxon_w_names_and_phenos_unique_combs=merge(x=mean_dev_by_taxon_w_names_unique_combs, y=phenos, by.x="AmesName", by.y="TaxaName") 
  rm(mean_dev_by_taxon_w_names_and_phenos_unique_combs)
  mean_dev_by_taxon_w_names_and_phenos_unique_combs=merge(x=mean_dev_by_taxon_w_names_unique_combs, y=phenos, by.x="RawPhenotypeNames", by.y="TaxaName") 
  write.table(cbind(TaxaNames=mean_dev_by_taxon_w_names_and_phenos_unique_combs$RawPhenotypeNames,Mean_dysreg=log10(mean_dev_by_taxon_w_names_and_phenos_unique_combs$Mean_Dev), FitnessBlup=mean_dev_by_taxon_w_names_and_phenos_unique_combs$TotalKernelWeight0607summer_blup), (paste(plot_out_dir,xlab,"vs_0607_seed_wt_BLUP.txt", sep="")), quote = F, row.names = F)
  write.table(mean_dev_by_taxon_w_names_and_phenos_unique_combs, paste(plot_out_dir,xlab,"vs_0607_seed_wt_BLUP_fulltable.txt", sep=""), quote = F, row.names = F)
  
  
  rm(linesummary)
  rm(line)
  #orig plot(mean_dev_by_taxon_w_names_and_phenos_unique_combs$Mean_Dev, mean_dev_by_taxon_w_names_and_phenos_unique_combs$PHT_BLUP, cex=1.5, col="deeppink4", pch=16, xlab=xlab, ylab=ylab)
  pdf(paste(plot_out_dir,xlab,"vs_0607_seed_wt_BLUP.pdf", sep=""),width=11,height=11,paper='special') 
  plot(log10(mean_dev_by_taxon_w_names_and_phenos_unique_combs$Mean_Dev), mean_dev_by_taxon_w_names_and_phenos_unique_combs$TotalKernelWeight0607summer_blup, cex=1.5, col="goldenrod3", pch=16, xlab=xlab, ylab="Seed Weight BLUPs")
  #orig line=lm(mean_dev_by_taxon_w_names_and_phenos_unique_combs$PHT_BLUP~mean_dev_by_taxon_w_names_and_phenos_unique_combs$Mean_Dev)
  line=lm(mean_dev_by_taxon_w_names_and_phenos_unique_combs$TotalKernelWeight0607summer_blup~log10(mean_dev_by_taxon_w_names_and_phenos_unique_combs$Mean_Dev))
  abline(line)
  linesummary=summary(line)
  #legend("topright", bty="n", legend=as.expression(c(bquote(R^2 == .(round(linesummary$r.squared, digits=4))), paste("P = ", bquote(.(round(pf(linesummary$fstatistic[1], linesummary$fstatistic[2], linesummary$fstatistic[3], lower.tail = FALSE), digits=10))), sep=""))))# p == .(pval1) .(pval2)))
  legend("topright", bty="n", legend=as.expression(c(bquote(r == .(round(sqrt(linesummary$r.squared), digits=4))), paste("P = ", bquote(.(round(pf(linesummary$fstatistic[1], linesummary$fstatistic[2], linesummary$fstatistic[3], lower.tail = FALSE), digits=10))), sep=""))))# p == .(pval1) .(pval2)))
  intercept=round(linesummary$coefficients[1,1], 2)
  slope=round(linesummary$coefficients[2,1], 2)
  sign=NULL
  if (slope>=0){sign="+"}
  legend("top", bty="n", legend=paste("y = ", intercept, sign, slope,"x", sep=""))
  dev.off()
  
  #PHTBLUP
  #Load in phenotypes 
  rm(phenos)
  phenos=read.table("/home/kak268/work_data/maize_phenotypes/Ames_BLUPs/Merged_Ames_traits", header=T)
  rm(mean_dev_by_taxon_w_names_and_phenos_unique_combs)
  mean_dev_by_taxon_w_names_and_phenos_unique_combs=merge(x=mean_dev_by_taxon_w_names_unique_combs, y=phenos, by.x="AmesName", by.y="V1") 
  #write.table(cbind(mean_dev_by_taxon_w_names_and_phenos_unique_combs$AmesName,log10(mean_dev_by_taxon_w_names_and_phenos_unique_combs$Mean_Dev), mean_dev_by_taxon_w_names_and_phenos_unique_combs$PHT_BLUP), (paste(plot_out_dir,xlab,"vs_mean_PHT_BLUP.pdf", sep="")), quote = F, row.names = F)
  write.table(cbind(TaxaNames=mean_dev_by_taxon_w_names_and_phenos_unique_combs$AmesName,Mean_dysreg=log10(mean_dev_by_taxon_w_names_and_phenos_unique_combs$Mean_Dev), FitnessBlup=mean_dev_by_taxon_w_names_and_phenos_unique_combs$PHT_BLUP), (paste(plot_out_dir,xlab,"vs_mean_PHT_BLUP.txt", sep="")), quote = F, row.names = F)
  write.table(mean_dev_by_taxon_w_names_and_phenos_unique_combs, paste(plot_out_dir,xlab,"vs_mean_PHT_BLUP_fulltable.txt", sep=""), quote = F, row.names = F)
  
  
  rm(linesummary)
  rm(line)
  pdf(paste(plot_out_dir,xlab, "vs_mean_PHT_BLUP.pdf", sep=""),width=11,height=11,paper='special') 
  plot(log10(mean_dev_by_taxon_w_names_and_phenos_unique_combs$Mean_Dev), mean_dev_by_taxon_w_names_and_phenos_unique_combs$PHT_BLUP, cex=1.5, col="cyan4", pch=16, xlab=xlab, ylab="Plant Height BLUP")
  line=lm(mean_dev_by_taxon_w_names_and_phenos_unique_combs$PHT_BLUP~log10(mean_dev_by_taxon_w_names_and_phenos_unique_combs$Mean_Dev))
  abline(line)
  linesummary=summary(line)
  #legend("topright", bty="n", legend=as.expression(c(bquote(R^2 == .(round(linesummary$r.squared, digits=4))), paste("P = ", bquote(.(round(pf(linesummary$fstatistic[1], linesummary$fstatistic[2], linesummary$fstatistic[3], lower.tail = FALSE), digits=5))), sep=""))))# p == .(pval1) .(pval2)))
  legend("topright", bty="n", legend=as.expression(c(bquote(r == .(round(sqrt(linesummary$r.squared), digits=4))), paste("P = ", bquote(.(round(pf(linesummary$fstatistic[1], linesummary$fstatistic[2], linesummary$fstatistic[3], lower.tail = FALSE), digits=5))), sep=""))))# p == .(pval1) .(pval2)))
  intercept=round(linesummary$coefficients[1,1], 2)
  slope=round(linesummary$coefficients[2,1], 2)
  sign=NULL
  if (slope>=0){sign="+"}
  legend("top", bty="n", legend=paste("y = ", intercept, sign, slope,"x", sep=""))
  dev.off()
  
  #GCA_FROM_YIELD_BLUPs
  #Load in phenotypes 
  # rm(phenos)
  # phenos=read.table("/home/kak268/work_data/maize_phenotypes/Yield_BLUPs_on_Testers/Ames_averaged_hyb_Yield_BLUPs_GCA.txt", header=T)
  # phenos=phenos[complete.cases(phenos),]
  # rm(mean_dev_by_taxon_w_names_and_phenos_unique_combs)
  # mean_dev_by_taxon_w_names_and_phenos_unique_combs=merge(x=mean_dev_by_taxon_w_names_unique_combs, y=phenos, by.x="AmesName", by.y="AmesName") 
  # mean_dev_by_taxon_w_names_and_phenos_unique_combs_by_RNAname=merge(x=mean_dev_by_taxon_w_names_unique_combs, y=phenos, by.x="RNA_TaxaName" , by.y="AmesName")
  # mean_dev_by_taxon_w_names_and_phenos_unique_combs_total=rbind(mean_dev_by_taxon_w_names_and_phenos_unique_combs_by_RNAname, mean_dev_by_taxon_w_names_and_phenos_unique_combs)
  # print(dim(mean_dev_by_taxon_w_names_and_phenos_unique_combs))
  # print(dim(mean_dev_by_taxon_w_names_and_phenos_unique_combs_by_RNAname))
  # print(dim(mean_dev_by_taxon_w_names_and_phenos_unique_combs_total))
  # #print(mean_dev_by_taxon_w_names_and_phenos_unique_combs_total)
  # 
  # rm(linesummary)
  # rm(line)
  # pdf(paste(plot_out_dir,xlab,"vs_mean_yield_BLUP.pdf", sep=""),width=11,height=11,paper='special') 
  # plot(log10(mean_dev_by_taxon_w_names_and_phenos_unique_combs_total$Mean_Dev), mean_dev_by_taxon_w_names_and_phenos_unique_combs_total$Mean_yield_BLUP, cex=1.5, col="darkred", pch=16, xlab=xlab, ylab="Mean Yield BLUP")
  # line=lm(mean_dev_by_taxon_w_names_and_phenos_unique_combs_total$Mean_yield_BLUP~log10(mean_dev_by_taxon_w_names_and_phenos_unique_combs_total$Mean_Dev))
  # abline(line)
  # linesummary=summary(line)
  # #legend("topright", bty="n", legend=as.expression(c(bquote(R^2 == .(round(linesummary$r.squared, digits=4))), paste("P = ", bquote(.(round(pf(linesummary$fstatistic[1], linesummary$fstatistic[2], linesummary$fstatistic[3], lower.tail = FALSE), digits=5))), sep=""))))# p == .(pval1) .(pval2)))
  # legend("topright", bty="n", legend=as.expression(c(bquote(r == .(round(sqrt(linesummary$r.squared), digits=4))), paste("P = ", bquote(.(round(pf(linesummary$fstatistic[1], linesummary$fstatistic[2], linesummary$fstatistic[3], lower.tail = FALSE), digits=5))), sep=""))))# p == .(pval1) .(pval2)))
  # intercept=round(linesummary$coefficients[1,1], 2)
  # slope=round(linesummary$coefficients[2,1], 2)
  # sign=NULL
  # if (slope>=0){sign="+"}
  # legend("top", bty="n", legend=paste("y = ", intercept, sign, slope,"x", sep=""))
  # dev.off()
}
#with median
dev_mat_master_median=NULL
plot_exp_dev_v_phenos_median=function(gene_w_exr_values, tissue, xlab, plot_out_dir){
  
  if(tissue!="all_tissues"){
    #Mean dev over all genes in each taxon
    #MEAN
    #mean_dev_by_taxon=cbind.data.frame(HMP32Name=as.character(var_mat$HMP32Name), Mean_Dev=as.numeric(rowMeans(var_mat[,2:ncol(var_mat)])))
    #MEDIAN
    #mean_dev_by_taxon=cbind.data.frame(HMP32Name=as.character(var_mat$HMP32Name), Mean_Dev=as.numeric(rowMedians(var_mat[,2:ncol(var_mat)])))
    #GEOMETRIC MEAN
    #mean_dev_by_taxon=cbind.data.frame(HMP32Name=as.character(var_mat$HMP32Name), Mean_Dev=as.numeric(apply(var_mat[,2:ncol(var_mat)], 1, function(x) exp(mean(log(x), na.rm=T)))))
    
    
    #Remove genes w 0 exp
    #To ensure that distantly related individuals which simply lack some B73 genes don't have inflated deviations, remove the variances from those genes which have zero expression
    var_mat_zero_exp_genes_set_to_zero_var=var_mat
    #var_mat_zero_exp_genes_set_to_zero_var[gene_w_exr_values==0]=NA
    
    #MEAN
    #mean_dev_by_taxon_exc_zero_exp_genes=cbind.data.frame(HMP32Name=as.character(var_mat_zero_exp_genes_set_to_zero_var$HMP32Name), Mean_Dev=as.numeric(rowMeans(var_mat_zero_exp_genes_set_to_zero_var[,2:ncol(var_mat_zero_exp_genes_set_to_zero_var)], na.rm=T)))
    #MEDIAN
    mean_dev_by_taxon_exc_zero_exp_genes=cbind.data.frame(HMP32Name=as.character(var_mat_zero_exp_genes_set_to_zero_var$HMP32Name), Mean_Dev=as.numeric(rowMedians(var_mat_zero_exp_genes_set_to_zero_var[,2:ncol(var_mat_zero_exp_genes_set_to_zero_var)], na.rm=T)))
    #GEOMETRIC MEAN
    #mean_dev_by_taxon_exc_zero_exp_genes=cbind.data.frame(HMP32Name=as.character(var_mat_zero_exp_genes_set_to_zero_var$HMP32Name), Mean_Dev=as.numeric(apply(var_mat_zero_exp_genes_set_to_zero_var[,2:ncol(var_mat_zero_exp_genes_set_to_zero_var)], 1, function(x) exp(mean(log(x), na.rm=T)))))
    
    #Seed Weight BLUPs
    #Load in phenotypes 
    #phenos=read.table("/home/kak268/work_data/maize_phenotypes/Ames_BLUPs/Merged_Ames_traits", header=T)
    phenos=read.table(file="/home/kak268/work_data/maize_phenotypes/282_panzea_phenotypes/traitMatrix_maize282NAM_v15-130212_Earweight_minus_Cobweight_melt_Means_and_BLUPS.txt", sep="\t", header=T)
    
    
    #Load in name converter so names canbe converted to Ames names
    name_decoder=read.table(header=T,"/media/kak268/A_2TB_Internal/RNA/Expressions_quants/samtools_idxstats_counts_bwamem_against_Zea_mays.AGPv3.29.cdna.all_w_Trinity_exc_NonAngiospermae_exc_over_85pct_B73_hits.fa/count_mat/renaming_scratch_folder/orig_1960names_and_alternate_names.txt")  
    name_decoder_tissue=name_decoder[name_decoder$TissueWODate==tissue,]
    ###Make the decoder unique
    unique_name_decoder_tissue=name_decoder_tissue[!duplicated(name_decoder_tissue$HMP32Name),]
    
    
    #Merge dev by taxon with name decoder and phenos
    #inc zero exp gene dev
    #mean_dev_by_taxon_w_names=merge(mean_dev_by_taxon, unique_name_decoder_tissue, by="HMP32Name")
    #exc zero exp gene dev
    mean_dev_by_taxon_w_names=merge(mean_dev_by_taxon_exc_zero_exp_genes, unique_name_decoder_tissue, by="HMP32Name")
    
    mean_dev_by_taxon_w_names_unique_combs=mean_dev_by_taxon_w_names[!duplicated(mean_dev_by_taxon_w_names[c("HMP32Name", "AmesName", "Mean_Dev", "RNA_TaxaName", "RawPhenotypeNames")]),]
    dev_mat_master_median=rbind(dev_mat_master_median, mean_dev_by_taxon_w_names_unique_combs)
    assign('dev_mat_master_median', dev_mat_master_median, envir=.GlobalEnv)
  } else{ # if tissue=="all_tissues"
    mean_dev_by_taxon_w_names_unique_combs=gene_w_exr_values    # in this case this holds the precomputed deviations from each tissue
  }
  
  #Seed Weight BLUPs
  #Load in phenotypes 
  #phenos=read.table("/home/kak268/work_data/maize_phenotypes/Ames_BLUPs/Merged_Ames_traits", header=T)
  phenos=read.table(file="/home/kak268/work_data/maize_phenotypes/282_panzea_phenotypes/traitMatrix_maize282NAM_v15-130212_Earweight_minus_Cobweight_melt_Means_and_BLUPS.txt", sep="\t", header=T)
  
  #orig#mean_dev_by_taxon_w_names_and_phenos_unique_combs=merge(x=mean_dev_by_taxon_w_names_unique_combs, y=phenos, by.x="AmesName", by.y="TaxaName") 
  mean_dev_by_taxon_w_names_and_phenos_unique_combs=merge(x=mean_dev_by_taxon_w_names_unique_combs, y=phenos, by.x="RawPhenotypeNames", by.y="TaxaName") 
  
  
  #orig plot(mean_dev_by_taxon_w_names_and_phenos_unique_combs$Mean_Dev, mean_dev_by_taxon_w_names_and_phenos_unique_combs$PHT_BLUP, cex=1.5, col="deeppink4", pch=16, xlab=xlab, ylab=ylab)
  pdf(paste(plot_out_dir,xlab,"vs_mean_seed_wt_BLUP.pdf", sep=""),width=11,height=11,paper='special') 
  plot(log10(mean_dev_by_taxon_w_names_and_phenos_unique_combs$Mean_Dev), mean_dev_by_taxon_w_names_and_phenos_unique_combs$SeedWeightBLUPs, cex=1.5, col="darkorange3", pch=16, xlab=xlab, ylab="Seed Weight BLUPs")
  #orig line=lm(mean_dev_by_taxon_w_names_and_phenos_unique_combs$PHT_BLUP~mean_dev_by_taxon_w_names_and_phenos_unique_combs$Mean_Dev)
  line=lm(mean_dev_by_taxon_w_names_and_phenos_unique_combs$SeedWeightBLUPs~log10(mean_dev_by_taxon_w_names_and_phenos_unique_combs$Mean_Dev))
  abline(line)
  linesummary=summary(line)
  #legend("topright", bty="n", legend=as.expression(c(bquote(R^2 == .(round(linesummary$r.squared, digits=4))), paste("P = ", bquote(.(round(pf(linesummary$fstatistic[1], linesummary$fstatistic[2], linesummary$fstatistic[3], lower.tail = FALSE), digits=10))), sep=""))))# p == .(pval1) .(pval2)))
  legend("topright", bty="n", legend=as.expression(c(bquote(r == .(round(sqrt(linesummary$r.squared), digits=4))), paste("P = ", bquote(.(round(pf(linesummary$fstatistic[1], linesummary$fstatistic[2], linesummary$fstatistic[3], lower.tail = FALSE), digits=10))), sep=""))))# p == .(pval1) .(pval2)))
  intercept=round(linesummary$coefficients[1,1], 2)
  slope=round(linesummary$coefficients[2,1], 2)
  sign=NULL
  if (slope>=0){sign="+"}
  legend("top", bty="n", legend=paste("y = ", intercept, sign, slope,"x", sep=""))
  dev.off()
  
  #HOLLAND 06 07 Seed Weight BLUPs
  rm(phenos)
  #Load in phenotypes 
  #phenos=read.table("/home/kak268/work_data/maize_phenotypes/Ames_BLUPs/Merged_Ames_traits", header=T)
  phenos=read.table(file="/home/kak268/work_data/maize_phenotypes/Holland_NAM_and_282_BLUPs/282 BLUPs 0607 all traits.txt", sep="\t", header=T)
  
  #orig#mean_dev_by_taxon_w_names_and_phenos_unique_combs=merge(x=mean_dev_by_taxon_w_names_unique_combs, y=phenos, by.x="AmesName", by.y="TaxaName") 
  rm(mean_dev_by_taxon_w_names_and_phenos_unique_combs)
  mean_dev_by_taxon_w_names_and_phenos_unique_combs=merge(x=mean_dev_by_taxon_w_names_unique_combs, y=phenos, by.x="RawPhenotypeNames", by.y="TaxaName") 
  
  rm(linesummary)
  rm(line)
  #orig plot(mean_dev_by_taxon_w_names_and_phenos_unique_combs$Mean_Dev, mean_dev_by_taxon_w_names_and_phenos_unique_combs$PHT_BLUP, cex=1.5, col="deeppink4", pch=16, xlab=xlab, ylab=ylab)
  pdf(paste(plot_out_dir,xlab,"vs_0607_seed_wt_BLUP.pdf", sep=""),width=11,height=11,paper='special') 
  plot(log10(mean_dev_by_taxon_w_names_and_phenos_unique_combs$Mean_Dev), mean_dev_by_taxon_w_names_and_phenos_unique_combs$TotalKernelWeight0607summer_blup, cex=1.5, col="goldenrod3", pch=16, xlab=xlab, ylab="Seed Weight BLUPs")
  #orig line=lm(mean_dev_by_taxon_w_names_and_phenos_unique_combs$PHT_BLUP~mean_dev_by_taxon_w_names_and_phenos_unique_combs$Mean_Dev)
  line=lm(mean_dev_by_taxon_w_names_and_phenos_unique_combs$TotalKernelWeight0607summer_blup~log10(mean_dev_by_taxon_w_names_and_phenos_unique_combs$Mean_Dev))
  abline(line)
  linesummary=summary(line)
  #legend("topright", bty="n", legend=as.expression(c(bquote(R^2 == .(round(linesummary$r.squared, digits=4))), paste("P = ", bquote(.(round(pf(linesummary$fstatistic[1], linesummary$fstatistic[2], linesummary$fstatistic[3], lower.tail = FALSE), digits=10))), sep=""))))# p == .(pval1) .(pval2)))
  legend("topright", bty="n", legend=as.expression(c(bquote(r == .(round(sqrt(linesummary$r.squared), digits=4))), paste("P = ", bquote(.(round(pf(linesummary$fstatistic[1], linesummary$fstatistic[2], linesummary$fstatistic[3], lower.tail = FALSE), digits=10))), sep=""))))# p == .(pval1) .(pval2)))
  intercept=round(linesummary$coefficients[1,1], 2)
  slope=round(linesummary$coefficients[2,1], 2)
  sign=NULL
  if (slope>=0){sign="+"}
  legend("top", bty="n", legend=paste("y = ", intercept, sign, slope,"x", sep=""))
  dev.off()
  
  #PHTBLUP
  #Load in phenotypes 
  rm(phenos)
  phenos=read.table("/home/kak268/work_data/maize_phenotypes/Ames_BLUPs/Merged_Ames_traits", header=T)
  rm(mean_dev_by_taxon_w_names_and_phenos_unique_combs)
  mean_dev_by_taxon_w_names_and_phenos_unique_combs=merge(x=mean_dev_by_taxon_w_names_unique_combs, y=phenos, by.x="AmesName", by.y="V1") 
  
  rm(linesummary)
  rm(line)
  pdf(paste(plot_out_dir,xlab, "vs_mean_PHT_BLUP.pdf", sep=""),width=11,height=11,paper='special') 
  plot(log10(mean_dev_by_taxon_w_names_and_phenos_unique_combs$Mean_Dev), mean_dev_by_taxon_w_names_and_phenos_unique_combs$PHT_BLUP, cex=1.5, col="cyan4", pch=16, xlab=xlab, ylab="Plant Height BLUP")
  line=lm(mean_dev_by_taxon_w_names_and_phenos_unique_combs$PHT_BLUP~log10(mean_dev_by_taxon_w_names_and_phenos_unique_combs$Mean_Dev))
  abline(line)
  linesummary=summary(line)
  #legend("topright", bty="n", legend=as.expression(c(bquote(R^2 == .(round(linesummary$r.squared, digits=4))), paste("P = ", bquote(.(round(pf(linesummary$fstatistic[1], linesummary$fstatistic[2], linesummary$fstatistic[3], lower.tail = FALSE), digits=5))), sep=""))))# p == .(pval1) .(pval2)))
  legend("topright", bty="n", legend=as.expression(c(bquote(r == .(round(sqrt(linesummary$r.squared), digits=4))), paste("P = ", bquote(.(round(pf(linesummary$fstatistic[1], linesummary$fstatistic[2], linesummary$fstatistic[3], lower.tail = FALSE), digits=5))), sep=""))))# p == .(pval1) .(pval2)))
  intercept=round(linesummary$coefficients[1,1], 2)
  slope=round(linesummary$coefficients[2,1], 2)
  sign=NULL
  if (slope>=0){sign="+"}
  legend("top", bty="n", legend=paste("y = ", intercept, sign, slope,"x", sep=""))
  dev.off()
  
  #GCA_FROM_YIELD_BLUPs
  #Load in phenotypes 
  rm(phenos)
  phenos=read.table("/home/kak268/work_data/maize_phenotypes/Yield_BLUPs_on_Testers/Ames_averaged_hyb_Yield_BLUPs_GCA.txt", header=T)
  phenos=phenos[complete.cases(phenos),]
  rm(mean_dev_by_taxon_w_names_and_phenos_unique_combs)
  mean_dev_by_taxon_w_names_and_phenos_unique_combs=merge(x=mean_dev_by_taxon_w_names_unique_combs, y=phenos, by.x="AmesName", by.y="AmesName") 
  mean_dev_by_taxon_w_names_and_phenos_unique_combs_by_RNAname=merge(x=mean_dev_by_taxon_w_names_unique_combs, y=phenos, by.x="RNA_TaxaName" , by.y="AmesName")
  mean_dev_by_taxon_w_names_and_phenos_unique_combs_total=rbind(mean_dev_by_taxon_w_names_and_phenos_unique_combs_by_RNAname, mean_dev_by_taxon_w_names_and_phenos_unique_combs)
  print(dim(mean_dev_by_taxon_w_names_and_phenos_unique_combs))
  print(dim(mean_dev_by_taxon_w_names_and_phenos_unique_combs_by_RNAname))
  print(dim(mean_dev_by_taxon_w_names_and_phenos_unique_combs_total))
  #print(mean_dev_by_taxon_w_names_and_phenos_unique_combs_total)
  
  rm(linesummary)
  rm(line)
  pdf(paste(plot_out_dir,xlab,"vs_mean_yield_BLUP.pdf", sep=""),width=11,height=11,paper='special') 
  plot(log10(mean_dev_by_taxon_w_names_and_phenos_unique_combs_total$Mean_Dev), mean_dev_by_taxon_w_names_and_phenos_unique_combs_total$Mean_yield_BLUP, cex=1.5, col="darkred", pch=16, xlab=xlab, ylab="Mean Yield BLUP")
  line=lm(mean_dev_by_taxon_w_names_and_phenos_unique_combs_total$Mean_yield_BLUP~log10(mean_dev_by_taxon_w_names_and_phenos_unique_combs_total$Mean_Dev))
  abline(line)
  linesummary=summary(line)
  #legend("topright", bty="n", legend=as.expression(c(bquote(R^2 == .(round(linesummary$r.squared, digits=4))), paste("P = ", bquote(.(round(pf(linesummary$fstatistic[1], linesummary$fstatistic[2], linesummary$fstatistic[3], lower.tail = FALSE), digits=5))), sep=""))))# p == .(pval1) .(pval2)))
  legend("topright", bty="n", legend=as.expression(c(bquote(r == .(round(sqrt(linesummary$r.squared), digits=4))), paste("P = ", bquote(.(round(pf(linesummary$fstatistic[1], linesummary$fstatistic[2], linesummary$fstatistic[3], lower.tail = FALSE), digits=5))), sep=""))))# p == .(pval1) .(pval2)))
  intercept=round(linesummary$coefficients[1,1], 2)
  slope=round(linesummary$coefficients[2,1], 2)
  sign=NULL
  if (slope>=0){sign="+"}
  legend("top", bty="n", legend=paste("y = ", intercept, sign, slope,"x", sep=""))
  dev.off()
}
#with geometric mean
dev_mat_master_geo_mean=NULL
plot_exp_dev_v_phenos_geo_mean=function(gene_w_exr_values, tissue, xlab, plot_out_dir){
  
  if(tissue!="all_tissues"){
    #Mean dev over all genes in each taxon
    #MEAN
    #mean_dev_by_taxon=cbind.data.frame(HMP32Name=as.character(var_mat$HMP32Name), Mean_Dev=as.numeric(rowMeans(var_mat[,2:ncol(var_mat)])))
    #MEDIAN
    #mean_dev_by_taxon=cbind.data.frame(HMP32Name=as.character(var_mat$HMP32Name), Mean_Dev=as.numeric(rowMedians(var_mat[,2:ncol(var_mat)])))
    #GEOMETRIC MEAN
    #mean_dev_by_taxon=cbind.data.frame(HMP32Name=as.character(var_mat$HMP32Name), Mean_Dev=as.numeric(apply(var_mat[,2:ncol(var_mat)], 1, function(x) exp(mean(log(x), na.rm=T)))))
    
    
    #Remove genes w 0 exp
    #To ensure that distantly related individuals which simply lack some B73 genes don't have inflated deviations, remove the variances from those genes which have zero expression
    var_mat_zero_exp_genes_set_to_zero_var=var_mat
    #var_mat_zero_exp_genes_set_to_zero_var[gene_w_exr_values==0]=NA
    
    #MEAN
    #mean_dev_by_taxon_exc_zero_exp_genes=cbind.data.frame(HMP32Name=as.character(var_mat_zero_exp_genes_set_to_zero_var$HMP32Name), Mean_Dev=as.numeric(rowMeans(var_mat_zero_exp_genes_set_to_zero_var[,2:ncol(var_mat_zero_exp_genes_set_to_zero_var)], na.rm=T)))
    #MEDIAN
    #mean_dev_by_taxon_exc_zero_exp_genes=cbind.data.frame(HMP32Name=as.character(var_mat_zero_exp_genes_set_to_zero_var$HMP32Name), Mean_Dev=as.numeric(rowMedians(var_mat_zero_exp_genes_set_to_zero_var[,2:ncol(var_mat_zero_exp_genes_set_to_zero_var)], na.rm=T)))
    #GEOMETRIC MEAN
    mean_dev_by_taxon_exc_zero_exp_genes=cbind.data.frame(HMP32Name=as.character(var_mat_zero_exp_genes_set_to_zero_var$HMP32Name), Mean_Dev=as.numeric(apply(var_mat_zero_exp_genes_set_to_zero_var[,2:ncol(var_mat_zero_exp_genes_set_to_zero_var)], 1, function(x) exp(mean(log(x), na.rm=T)))))
    
    #Seed Weight BLUPs
    #Load in phenotypes 
    #phenos=read.table("/home/kak268/work_data/maize_phenotypes/Ames_BLUPs/Merged_Ames_traits", header=T)
    phenos=read.table(file="/home/kak268/work_data/maize_phenotypes/282_panzea_phenotypes/traitMatrix_maize282NAM_v15-130212_Earweight_minus_Cobweight_melt_Means_and_BLUPS.txt", sep="\t", header=T)
    
    
    #Load in name converter so names canbe converted to Ames names
    name_decoder=read.table(header=T,"/media/kak268/A_2TB_Internal/RNA/Expressions_quants/samtools_idxstats_counts_bwamem_against_Zea_mays.AGPv3.29.cdna.all_w_Trinity_exc_NonAngiospermae_exc_over_85pct_B73_hits.fa/count_mat/renaming_scratch_folder/orig_1960names_and_alternate_names.txt")  
    name_decoder_tissue=name_decoder[name_decoder$TissueWODate==tissue,]
    ###Make the decoder unique
    unique_name_decoder_tissue=name_decoder_tissue[!duplicated(name_decoder_tissue$HMP32Name),]
    
    
    #Merge dev by taxon with name decoder and phenos
    #inc zero exp gene dev
    #mean_dev_by_taxon_w_names=merge(mean_dev_by_taxon, unique_name_decoder_tissue, by="HMP32Name")
    #exc zero exp gene dev
    mean_dev_by_taxon_w_names=merge(mean_dev_by_taxon_exc_zero_exp_genes, unique_name_decoder_tissue, by="HMP32Name")
    
    mean_dev_by_taxon_w_names_unique_combs=mean_dev_by_taxon_w_names[!duplicated(mean_dev_by_taxon_w_names[c("HMP32Name", "AmesName", "Mean_Dev", "RNA_TaxaName", "RawPhenotypeNames")]),]
    dev_mat_master_geo_mean=rbind(dev_mat_master_geo_mean, mean_dev_by_taxon_w_names_unique_combs)
    assign('dev_mat_master_geo_mean', dev_mat_master_geo_mean, envir=.GlobalEnv)
  } else{ # if tissue=="all_tissues"
    mean_dev_by_taxon_w_names_unique_combs=gene_w_exr_values    # in this case this holds the precomputed deviations from each tissue
  }
  
  #Seed Weight BLUPs
  #Load in phenotypes 
  #phenos=read.table("/home/kak268/work_data/maize_phenotypes/Ames_BLUPs/Merged_Ames_traits", header=T)
  phenos=read.table(file="/home/kak268/work_data/maize_phenotypes/282_panzea_phenotypes/traitMatrix_maize282NAM_v15-130212_Earweight_minus_Cobweight_melt_Means_and_BLUPS.txt", sep="\t", header=T)
  
  #orig#mean_dev_by_taxon_w_names_and_phenos_unique_combs=merge(x=mean_dev_by_taxon_w_names_unique_combs, y=phenos, by.x="AmesName", by.y="TaxaName") 
  mean_dev_by_taxon_w_names_and_phenos_unique_combs=merge(x=mean_dev_by_taxon_w_names_unique_combs, y=phenos, by.x="RawPhenotypeNames", by.y="TaxaName") 
  
  
  #orig plot(mean_dev_by_taxon_w_names_and_phenos_unique_combs$Mean_Dev, mean_dev_by_taxon_w_names_and_phenos_unique_combs$PHT_BLUP, cex=1.5, col="deeppink4", pch=16, xlab=xlab, ylab=ylab)
  pdf(paste(plot_out_dir,xlab,"vs_mean_seed_wt_BLUP.pdf", sep=""),width=11,height=11,paper='special') 
  plot(log10(mean_dev_by_taxon_w_names_and_phenos_unique_combs$Mean_Dev), mean_dev_by_taxon_w_names_and_phenos_unique_combs$SeedWeightBLUPs, cex=1.5, col="darkorange3", pch=16, xlab=xlab, ylab="Seed Weight BLUPs")
  #orig line=lm(mean_dev_by_taxon_w_names_and_phenos_unique_combs$PHT_BLUP~mean_dev_by_taxon_w_names_and_phenos_unique_combs$Mean_Dev)
  line=lm(mean_dev_by_taxon_w_names_and_phenos_unique_combs$SeedWeightBLUPs~log10(mean_dev_by_taxon_w_names_and_phenos_unique_combs$Mean_Dev))
  abline(line)
  linesummary=summary(line)
  #legend("topright", bty="n", legend=as.expression(c(bquote(R^2 == .(round(linesummary$r.squared, digits=4))), paste("P = ", bquote(.(round(pf(linesummary$fstatistic[1], linesummary$fstatistic[2], linesummary$fstatistic[3], lower.tail = FALSE), digits=10))), sep=""))))# p == .(pval1) .(pval2)))
  legend("topright", bty="n", legend=as.expression(c(bquote(r == .(round(sqrt(linesummary$r.squared), digits=4))), paste("P = ", bquote(.(round(pf(linesummary$fstatistic[1], linesummary$fstatistic[2], linesummary$fstatistic[3], lower.tail = FALSE), digits=10))), sep=""))))# p == .(pval1) .(pval2)))
  intercept=round(linesummary$coefficients[1,1], 2)
  slope=round(linesummary$coefficients[2,1], 2)
  sign=NULL
  if (slope>=0){sign="+"}
  legend("top", bty="n", legend=paste("y = ", intercept, sign, slope,"x", sep=""))
  dev.off()
  
  #HOLLAND 06 07 Seed Weight BLUPs
  rm(phenos)
  #Load in phenotypes 
  #phenos=read.table("/home/kak268/work_data/maize_phenotypes/Ames_BLUPs/Merged_Ames_traits", header=T)
  phenos=read.table(file="/home/kak268/work_data/maize_phenotypes/Holland_NAM_and_282_BLUPs/282 BLUPs 0607 all traits.txt", sep="\t", header=T)
  
  #orig#mean_dev_by_taxon_w_names_and_phenos_unique_combs=merge(x=mean_dev_by_taxon_w_names_unique_combs, y=phenos, by.x="AmesName", by.y="TaxaName") 
  rm(mean_dev_by_taxon_w_names_and_phenos_unique_combs)
  mean_dev_by_taxon_w_names_and_phenos_unique_combs=merge(x=mean_dev_by_taxon_w_names_unique_combs, y=phenos, by.x="RawPhenotypeNames", by.y="TaxaName") 
  
  rm(linesummary)
  rm(line)
  #orig plot(mean_dev_by_taxon_w_names_and_phenos_unique_combs$Mean_Dev, mean_dev_by_taxon_w_names_and_phenos_unique_combs$PHT_BLUP, cex=1.5, col="deeppink4", pch=16, xlab=xlab, ylab=ylab)
  pdf(paste(plot_out_dir,xlab,"vs_0607_seed_wt_BLUP.pdf", sep=""),width=11,height=11,paper='special') 
  plot(log10(mean_dev_by_taxon_w_names_and_phenos_unique_combs$Mean_Dev), mean_dev_by_taxon_w_names_and_phenos_unique_combs$TotalKernelWeight0607summer_blup, cex=1.5, col="goldenrod3", pch=16, xlab=xlab, ylab="Seed Weight BLUPs")
  #orig line=lm(mean_dev_by_taxon_w_names_and_phenos_unique_combs$PHT_BLUP~mean_dev_by_taxon_w_names_and_phenos_unique_combs$Mean_Dev)
  line=lm(mean_dev_by_taxon_w_names_and_phenos_unique_combs$TotalKernelWeight0607summer_blup~log10(mean_dev_by_taxon_w_names_and_phenos_unique_combs$Mean_Dev))
  abline(line)
  linesummary=summary(line)
  #legend("topright", bty="n", legend=as.expression(c(bquote(R^2 == .(round(linesummary$r.squared, digits=4))), paste("P = ", bquote(.(round(pf(linesummary$fstatistic[1], linesummary$fstatistic[2], linesummary$fstatistic[3], lower.tail = FALSE), digits=10))), sep=""))))# p == .(pval1) .(pval2)))
  legend("topright", bty="n", legend=as.expression(c(bquote(r == .(round(sqrt(linesummary$r.squared), digits=4))), paste("P = ", bquote(.(round(pf(linesummary$fstatistic[1], linesummary$fstatistic[2], linesummary$fstatistic[3], lower.tail = FALSE), digits=10))), sep=""))))# p == .(pval1) .(pval2)))
  intercept=round(linesummary$coefficients[1,1], 2)
  slope=round(linesummary$coefficients[2,1], 2)
  sign=NULL
  if (slope>=0){sign="+"}
  legend("top", bty="n", legend=paste("y = ", intercept, sign, slope,"x", sep=""))
  dev.off()
  
  #PHTBLUP
  #Load in phenotypes 
  rm(phenos)
  phenos=read.table("/home/kak268/work_data/maize_phenotypes/Ames_BLUPs/Merged_Ames_traits", header=T)
  rm(mean_dev_by_taxon_w_names_and_phenos_unique_combs)
  mean_dev_by_taxon_w_names_and_phenos_unique_combs=merge(x=mean_dev_by_taxon_w_names_unique_combs, y=phenos, by.x="AmesName", by.y="V1") 
  
  rm(linesummary)
  rm(line)
  pdf(paste(plot_out_dir,xlab, "vs_mean_PHT_BLUP.pdf", sep=""),width=11,height=11,paper='special') 
  plot(log10(mean_dev_by_taxon_w_names_and_phenos_unique_combs$Mean_Dev), mean_dev_by_taxon_w_names_and_phenos_unique_combs$PHT_BLUP, cex=1.5, col="cyan4", pch=16, xlab=xlab, ylab="Plant Height BLUP")
  line=lm(mean_dev_by_taxon_w_names_and_phenos_unique_combs$PHT_BLUP~log10(mean_dev_by_taxon_w_names_and_phenos_unique_combs$Mean_Dev))
  abline(line)
  linesummary=summary(line)
  #legend("topright", bty="n", legend=as.expression(c(bquote(R^2 == .(round(linesummary$r.squared, digits=4))), paste("P = ", bquote(.(round(pf(linesummary$fstatistic[1], linesummary$fstatistic[2], linesummary$fstatistic[3], lower.tail = FALSE), digits=5))), sep=""))))# p == .(pval1) .(pval2)))
  legend("topright", bty="n", legend=as.expression(c(bquote(r == .(round(sqrt(linesummary$r.squared), digits=4))), paste("P = ", bquote(.(round(pf(linesummary$fstatistic[1], linesummary$fstatistic[2], linesummary$fstatistic[3], lower.tail = FALSE), digits=5))), sep=""))))# p == .(pval1) .(pval2)))
  intercept=round(linesummary$coefficients[1,1], 2)
  slope=round(linesummary$coefficients[2,1], 2)
  sign=NULL
  if (slope>=0){sign="+"}
  legend("top", bty="n", legend=paste("y = ", intercept, sign, slope,"x", sep=""))
  dev.off()
  
  #GCA_FROM_YIELD_BLUPs
  #Load in phenotypes 
  rm(phenos)
  phenos=read.table("/home/kak268/work_data/maize_phenotypes/Yield_BLUPs_on_Testers/Ames_averaged_hyb_Yield_BLUPs_GCA.txt", header=T)
  phenos=phenos[complete.cases(phenos),]
  rm(mean_dev_by_taxon_w_names_and_phenos_unique_combs)
  mean_dev_by_taxon_w_names_and_phenos_unique_combs=merge(x=mean_dev_by_taxon_w_names_unique_combs, y=phenos, by.x="AmesName", by.y="AmesName") 
  mean_dev_by_taxon_w_names_and_phenos_unique_combs_by_RNAname=merge(x=mean_dev_by_taxon_w_names_unique_combs, y=phenos, by.x="RNA_TaxaName" , by.y="AmesName")
  mean_dev_by_taxon_w_names_and_phenos_unique_combs_total=rbind(mean_dev_by_taxon_w_names_and_phenos_unique_combs_by_RNAname, mean_dev_by_taxon_w_names_and_phenos_unique_combs)
  print(dim(mean_dev_by_taxon_w_names_and_phenos_unique_combs))
  print(dim(mean_dev_by_taxon_w_names_and_phenos_unique_combs_by_RNAname))
  print(dim(mean_dev_by_taxon_w_names_and_phenos_unique_combs_total))
  #print(mean_dev_by_taxon_w_names_and_phenos_unique_combs_total)
  
  rm(linesummary)
  rm(line)
  pdf(paste(plot_out_dir,xlab,"vs_mean_yield_BLUP.pdf", sep=""),width=11,height=11,paper='special') 
  plot(log10(mean_dev_by_taxon_w_names_and_phenos_unique_combs_total$Mean_Dev), mean_dev_by_taxon_w_names_and_phenos_unique_combs_total$Mean_yield_BLUP, cex=1.5, col="darkred", pch=16, xlab=xlab, ylab="Mean Yield BLUP")
  line=lm(mean_dev_by_taxon_w_names_and_phenos_unique_combs_total$Mean_yield_BLUP~log10(mean_dev_by_taxon_w_names_and_phenos_unique_combs_total$Mean_Dev))
  abline(line)
  linesummary=summary(line)
  #legend("topright", bty="n", legend=as.expression(c(bquote(R^2 == .(round(linesummary$r.squared, digits=4))), paste("P = ", bquote(.(round(pf(linesummary$fstatistic[1], linesummary$fstatistic[2], linesummary$fstatistic[3], lower.tail = FALSE), digits=5))), sep=""))))# p == .(pval1) .(pval2)))
  legend("topright", bty="n", legend=as.expression(c(bquote(r == .(round(sqrt(linesummary$r.squared), digits=4))), paste("P = ", bquote(.(round(pf(linesummary$fstatistic[1], linesummary$fstatistic[2], linesummary$fstatistic[3], lower.tail = FALSE), digits=5))), sep=""))))# p == .(pval1) .(pval2)))
  intercept=round(linesummary$coefficients[1,1], 2)
  slope=round(linesummary$coefficients[2,1], 2)
  sign=NULL
  if (slope>=0){sign="+"}
  legend("top", bty="n", legend=paste("y = ", intercept, sign, slope,"x", sep=""))
  dev.off()
}

##
#Rank based
dev_mat_rank_master_mean=NULL
plot_exp_dev_v_phenos_rank=function(gene_w_exr_values, tissue, xlab, plot_out_dir){
  
  if(tissue!="all_tissues"){
    #Mean dev over all genes in each taxon
    #MEAN
    #mean_dev_by_taxon=cbind.data.frame(HMP32Name=as.character(var_mat$HMP32Name), Mean_Dev=as.numeric(rowMeans(var_mat[,2:ncol(var_mat)])))
    #MEDIAN
    #mean_dev_by_taxon=cbind.data.frame(HMP32Name=as.character(var_mat$HMP32Name), Mean_Dev=as.numeric(rowMedians(var_mat[,2:ncol(var_mat)])))
    #GEOMETRIC MEAN
    #mean_dev_by_taxon=cbind.data.frame(HMP32Name=as.character(var_mat$HMP32Name), Mean_Dev=as.numeric(apply(var_mat[,2:ncol(var_mat)], 1, function(x) exp(mean(log(x), na.rm=T)))))
    
    
    #Remove genes w 0 exp
    #To ensure that distantly related individuals which simply lack some B73 genes don't have inflated deviations, remove the variances from those genes which have zero expression
    var_mat_zero_exp_genes_set_to_zero_var=var_mat_rank
    #var_mat_zero_exp_genes_set_to_zero_var[gene_w_exr_values==0]=NA
    
    #MEAN
    mean_dev_by_taxon_exc_zero_exp_genes=cbind.data.frame(HMP32Name=as.character(var_mat_zero_exp_genes_set_to_zero_var$HMP32Name), Mean_Dev=as.numeric(rowMeans(var_mat_zero_exp_genes_set_to_zero_var[,2:ncol(var_mat_zero_exp_genes_set_to_zero_var)], na.rm=T)))
    #MEDIAN
    #mean_dev_by_taxon_exc_zero_exp_genes=cbind.data.frame(HMP32Name=as.character(var_mat_zero_exp_genes_set_to_zero_var$HMP32Name), Mean_Dev=as.numeric(rowMedians(var_mat_zero_exp_genes_set_to_zero_var[,2:ncol(var_mat_zero_exp_genes_set_to_zero_var)], na.rm=T)))
    #GEOMETRIC MEAN
    #mean_dev_by_taxon_exc_zero_exp_genes=cbind.data.frame(HMP32Name=as.character(var_mat_zero_exp_genes_set_to_zero_var$HMP32Name), Mean_Dev=as.numeric(apply(var_mat_zero_exp_genes_set_to_zero_var[,2:ncol(var_mat_zero_exp_genes_set_to_zero_var)], 1, function(x) exp(mean(log(x), na.rm=T)))))
    
    #Seed Weight BLUPs
    #Load in phenotypes 
    #phenos=read.table("/home/kak268/work_data/maize_phenotypes/Ames_BLUPs/Merged_Ames_traits", header=T)
    phenos=read.table(file="/home/kak268/work_data/maize_phenotypes/282_panzea_phenotypes/traitMatrix_maize282NAM_v15-130212_Earweight_minus_Cobweight_melt_Means_and_BLUPS.txt", sep="\t", header=T)
    
    
    #Load in name converter so names canbe converted to Ames names
    name_decoder=read.table(header=T,"/media/kak268/A_2TB_Internal/RNA/Expressions_quants/samtools_idxstats_counts_bwamem_against_Zea_mays.AGPv3.29.cdna.all_w_Trinity_exc_NonAngiospermae_exc_over_85pct_B73_hits.fa/count_mat/renaming_scratch_folder/orig_1960names_and_alternate_names.txt")  
    name_decoder_tissue=name_decoder[name_decoder$TissueWODate==tissue,]
    ###Make the decoder unique
    unique_name_decoder_tissue=name_decoder_tissue[!duplicated(name_decoder_tissue$HMP32Name),]
    
    
    #Merge dev by taxon with name decoder and phenos
    #inc zero exp gene dev
    #mean_dev_by_taxon_w_names=merge(mean_dev_by_taxon, unique_name_decoder_tissue, by="HMP32Name")
    #exc zero exp gene dev
    mean_dev_by_taxon_w_names=merge(mean_dev_by_taxon_exc_zero_exp_genes, unique_name_decoder_tissue, by="HMP32Name")
    
    mean_dev_by_taxon_w_names_unique_combs=mean_dev_by_taxon_w_names[!duplicated(mean_dev_by_taxon_w_names[c("HMP32Name", "AmesName", "Mean_Dev", "RNA_TaxaName", "RawPhenotypeNames")]),]
    dev_mat_rank_master_mean=rbind(dev_mat_rank_master_mean, mean_dev_by_taxon_w_names_unique_combs)
    assign('dev_mat_rank_master_mean', dev_mat_rank_master_mean, envir=.GlobalEnv)
  } else{ # if tissue=="all_tissues"
    mean_dev_by_taxon_w_names_unique_combs=gene_w_exr_values    # in this case this holds the precomputed deviations from each tissue
  }
  
  #Seed Weight BLUPs
  #Load in phenotypes 
  #phenos=read.table("/home/kak268/work_data/maize_phenotypes/Ames_BLUPs/Merged_Ames_traits", header=T)
  phenos=read.table(file="/home/kak268/work_data/maize_phenotypes/282_panzea_phenotypes/traitMatrix_maize282NAM_v15-130212_Earweight_minus_Cobweight_melt_Means_and_BLUPS.txt", sep="\t", header=T)
  
  #orig#mean_dev_by_taxon_w_names_and_phenos_unique_combs=merge(x=mean_dev_by_taxon_w_names_unique_combs, y=phenos, by.x="AmesName", by.y="TaxaName") 
  mean_dev_by_taxon_w_names_and_phenos_unique_combs=merge(x=mean_dev_by_taxon_w_names_unique_combs, y=phenos, by.x="RawPhenotypeNames", by.y="TaxaName") 
  
  
  #orig plot(mean_dev_by_taxon_w_names_and_phenos_unique_combs$Mean_Dev, mean_dev_by_taxon_w_names_and_phenos_unique_combs$PHT_BLUP, cex=1.5, col="deeppink4", pch=16, xlab=xlab, ylab=ylab)
  pdf(paste(plot_out_dir,xlab,"vs_mean_seed_wt_BLUP.pdf", sep=""),width=11,height=11,paper='special') 
  plot(log10(mean_dev_by_taxon_w_names_and_phenos_unique_combs$Mean_Dev), mean_dev_by_taxon_w_names_and_phenos_unique_combs$SeedWeightBLUPs, cex=1.5, col="darkorange3", pch=16, xlab=xlab, ylab="Seed Weight BLUPs")
  #orig line=lm(mean_dev_by_taxon_w_names_and_phenos_unique_combs$PHT_BLUP~mean_dev_by_taxon_w_names_and_phenos_unique_combs$Mean_Dev)
  line=lm(mean_dev_by_taxon_w_names_and_phenos_unique_combs$SeedWeightBLUPs~log10(mean_dev_by_taxon_w_names_and_phenos_unique_combs$Mean_Dev))
  abline(line)
  linesummary=summary(line)
  #legend("topright", bty="n", legend=as.expression(c(bquote(R^2 == .(round(linesummary$r.squared, digits=4))), paste("P = ", bquote(.(round(pf(linesummary$fstatistic[1], linesummary$fstatistic[2], linesummary$fstatistic[3], lower.tail = FALSE), digits=10))), sep=""))))# p == .(pval1) .(pval2)))
  legend("topright", bty="n", legend=as.expression(c(bquote(r == .(round(sqrt(linesummary$r.squared), digits=4))), paste("P = ", bquote(.(round(pf(linesummary$fstatistic[1], linesummary$fstatistic[2], linesummary$fstatistic[3], lower.tail = FALSE), digits=10))), sep=""))))# p == .(pval1) .(pval2)))
  intercept=round(linesummary$coefficients[1,1], 2)
  slope=round(linesummary$coefficients[2,1], 2)
  sign=NULL
  if (slope>=0){sign="+"}
  legend("top", bty="n", legend=paste("y = ", intercept, sign, slope,"x", sep=""))
  dev.off()
  
  #HOLLAND 06 07 Seed Weight BLUPs
  rm(phenos)
  #Load in phenotypes 
  #phenos=read.table("/home/kak268/work_data/maize_phenotypes/Ames_BLUPs/Merged_Ames_traits", header=T)
  phenos=read.table(file="/home/kak268/work_data/maize_phenotypes/Holland_NAM_and_282_BLUPs/282 BLUPs 0607 all traits.txt", sep="\t", header=T)
  
  #orig#mean_dev_by_taxon_w_names_and_phenos_unique_combs=merge(x=mean_dev_by_taxon_w_names_unique_combs, y=phenos, by.x="AmesName", by.y="TaxaName") 
  rm(mean_dev_by_taxon_w_names_and_phenos_unique_combs)
  mean_dev_by_taxon_w_names_and_phenos_unique_combs=merge(x=mean_dev_by_taxon_w_names_unique_combs, y=phenos, by.x="RawPhenotypeNames", by.y="TaxaName") 
  
  rm(linesummary)
  rm(line)
  #orig plot(mean_dev_by_taxon_w_names_and_phenos_unique_combs$Mean_Dev, mean_dev_by_taxon_w_names_and_phenos_unique_combs$PHT_BLUP, cex=1.5, col="deeppink4", pch=16, xlab=xlab, ylab=ylab)
  pdf(paste(plot_out_dir,xlab,"vs_0607_seed_wt_BLUP.pdf", sep=""),width=11,height=11,paper='special') 
  plot(log10(mean_dev_by_taxon_w_names_and_phenos_unique_combs$Mean_Dev), mean_dev_by_taxon_w_names_and_phenos_unique_combs$TotalKernelWeight0607summer_blup, cex=1.5, col="goldenrod3", pch=16, xlab=xlab, ylab="Seed Weight BLUPs")
  #orig line=lm(mean_dev_by_taxon_w_names_and_phenos_unique_combs$PHT_BLUP~mean_dev_by_taxon_w_names_and_phenos_unique_combs$Mean_Dev)
  line=lm(mean_dev_by_taxon_w_names_and_phenos_unique_combs$TotalKernelWeight0607summer_blup~log10(mean_dev_by_taxon_w_names_and_phenos_unique_combs$Mean_Dev))
  abline(line)
  linesummary=summary(line)
  #legend("topright", bty="n", legend=as.expression(c(bquote(R^2 == .(round(linesummary$r.squared, digits=4))), paste("P = ", bquote(.(round(pf(linesummary$fstatistic[1], linesummary$fstatistic[2], linesummary$fstatistic[3], lower.tail = FALSE), digits=10))), sep=""))))# p == .(pval1) .(pval2)))
  legend("topright", bty="n", legend=as.expression(c(bquote(r == .(round(sqrt(linesummary$r.squared), digits=4))), paste("P = ", bquote(.(round(pf(linesummary$fstatistic[1], linesummary$fstatistic[2], linesummary$fstatistic[3], lower.tail = FALSE), digits=10))), sep=""))))# p == .(pval1) .(pval2)))
  intercept=round(linesummary$coefficients[1,1], 2)
  slope=round(linesummary$coefficients[2,1], 2)
  sign=NULL
  if (slope>=0){sign="+"}
  legend("top", bty="n", legend=paste("y = ", intercept, sign, slope,"x", sep=""))
  dev.off()
  
  #PHTBLUP
  #Load in phenotypes 
  rm(phenos)
  phenos=read.table("/home/kak268/work_data/maize_phenotypes/Ames_BLUPs/Merged_Ames_traits", header=T)
  rm(mean_dev_by_taxon_w_names_and_phenos_unique_combs)
  mean_dev_by_taxon_w_names_and_phenos_unique_combs=merge(x=mean_dev_by_taxon_w_names_unique_combs, y=phenos, by.x="AmesName", by.y="V1") 
  
  rm(linesummary)
  rm(line)
  pdf(paste(plot_out_dir,xlab, "vs_mean_PHT_BLUP.pdf", sep=""),width=11,height=11,paper='special') 
  plot(log10(mean_dev_by_taxon_w_names_and_phenos_unique_combs$Mean_Dev), mean_dev_by_taxon_w_names_and_phenos_unique_combs$PHT_BLUP, cex=1.5, col="cyan4", pch=16, xlab=xlab, ylab="Plant Height BLUP")
  line=lm(mean_dev_by_taxon_w_names_and_phenos_unique_combs$PHT_BLUP~log10(mean_dev_by_taxon_w_names_and_phenos_unique_combs$Mean_Dev))
  abline(line)
  linesummary=summary(line)
  #legend("topright", bty="n", legend=as.expression(c(bquote(R^2 == .(round(linesummary$r.squared, digits=4))), paste("P = ", bquote(.(round(pf(linesummary$fstatistic[1], linesummary$fstatistic[2], linesummary$fstatistic[3], lower.tail = FALSE), digits=5))), sep=""))))# p == .(pval1) .(pval2)))
  legend("topright", bty="n", legend=as.expression(c(bquote(r == .(round(sqrt(linesummary$r.squared), digits=4))), paste("P = ", bquote(.(round(pf(linesummary$fstatistic[1], linesummary$fstatistic[2], linesummary$fstatistic[3], lower.tail = FALSE), digits=5))), sep=""))))# p == .(pval1) .(pval2)))
  intercept=round(linesummary$coefficients[1,1], 2)
  slope=round(linesummary$coefficients[2,1], 2)
  sign=NULL
  if (slope>=0){sign="+"}
  legend("top", bty="n", legend=paste("y = ", intercept, sign, slope,"x", sep=""))
  dev.off()
  
  #GCA_FROM_YIELD_BLUPs
  #Load in phenotypes 
  rm(phenos)
  phenos=read.table("/home/kak268/work_data/maize_phenotypes/Yield_BLUPs_on_Testers/Ames_averaged_hyb_Yield_BLUPs_GCA.txt", header=T)
  phenos=phenos[complete.cases(phenos),]
  rm(mean_dev_by_taxon_w_names_and_phenos_unique_combs)
  mean_dev_by_taxon_w_names_and_phenos_unique_combs=merge(x=mean_dev_by_taxon_w_names_unique_combs, y=phenos, by.x="AmesName", by.y="AmesName") 
  mean_dev_by_taxon_w_names_and_phenos_unique_combs_by_RNAname=merge(x=mean_dev_by_taxon_w_names_unique_combs, y=phenos, by.x="RNA_TaxaName" , by.y="AmesName")
  mean_dev_by_taxon_w_names_and_phenos_unique_combs_total=rbind(mean_dev_by_taxon_w_names_and_phenos_unique_combs_by_RNAname, mean_dev_by_taxon_w_names_and_phenos_unique_combs)
  print(dim(mean_dev_by_taxon_w_names_and_phenos_unique_combs))
  print(dim(mean_dev_by_taxon_w_names_and_phenos_unique_combs_by_RNAname))
  print(dim(mean_dev_by_taxon_w_names_and_phenos_unique_combs_total))
  #print(mean_dev_by_taxon_w_names_and_phenos_unique_combs_total)
  
  rm(linesummary)
  rm(line)
  pdf(paste(plot_out_dir,xlab,"vs_mean_yield_BLUP.pdf", sep=""),width=11,height=11,paper='special') 
  plot(log10(mean_dev_by_taxon_w_names_and_phenos_unique_combs_total$Mean_Dev), mean_dev_by_taxon_w_names_and_phenos_unique_combs_total$Mean_yield_BLUP, cex=1.5, col="darkred", pch=16, xlab=xlab, ylab="Mean Yield BLUP")
  line=lm(mean_dev_by_taxon_w_names_and_phenos_unique_combs_total$Mean_yield_BLUP~log10(mean_dev_by_taxon_w_names_and_phenos_unique_combs_total$Mean_Dev))
  abline(line)
  linesummary=summary(line)
  #legend("topright", bty="n", legend=as.expression(c(bquote(R^2 == .(round(linesummary$r.squared, digits=4))), paste("P = ", bquote(.(round(pf(linesummary$fstatistic[1], linesummary$fstatistic[2], linesummary$fstatistic[3], lower.tail = FALSE), digits=5))), sep=""))))# p == .(pval1) .(pval2)))
  legend("topright", bty="n", legend=as.expression(c(bquote(r == .(round(sqrt(linesummary$r.squared), digits=4))), paste("P = ", bquote(.(round(pf(linesummary$fstatistic[1], linesummary$fstatistic[2], linesummary$fstatistic[3], lower.tail = FALSE), digits=5))), sep=""))))# p == .(pval1) .(pval2)))
  intercept=round(linesummary$coefficients[1,1], 2)
  slope=round(linesummary$coefficients[2,1], 2)
  sign=NULL
  if (slope>=0){sign="+"}
  legend("top", bty="n", legend=paste("y = ", intercept, sign, slope,"x", sep=""))
  dev.off()
}



plot_exp_dev_by_gene_v_phenos=function(gene_w_exr_values, xlab, plot_out_dir){
#   means=colMeans(gene_w_exr_values[,2:ncol(gene_w_exr_values)])
#   var_mat=gene_w_exr_values
#   var_mat[,-1]=NA
#   
#   for (i in 1:nrow(gene_w_exr_values)){
#     var_mat[i,2:ncol(var_mat)]=(gene_w_exr_values[i,2:ncol(var_mat)]-means)^2
#     cat(i)
#   }
    
  #Mean dev over all genes in each taxon
  #MEAN
  #mean_dev_by_taxon=cbind.data.frame(HMP32Name=as.character(var_mat$HMP32Name), Mean_Dev=as.numeric(rowMeans(var_mat[,2:ncol(var_mat)])))
  #MEDIAN
  #mean_dev_by_taxon=cbind.data.frame(HMP32Name=as.character(var_mat$HMP32Name), Mean_Dev=as.numeric(rowMedians(var_mat[,2:ncol(var_mat)])))
  #GEOMETRIC MEAN
  #mean_dev_by_taxon=cbind.data.frame(HMP32Name=as.character(var_mat$HMP32Name), Mean_Dev=as.numeric(apply(var_mat[,2:ncol(var_mat)], 1, function(x) exp(mean(log(x), na.rm=T)))))
  
  
  #Remove genes w 0 exp
  #To ensure that distantly related individuals which simply lack some B73 genes don't have inflated deviations, remove the variances from those genes which have zero expression
  var_mat_zero_exp_genes_set_to_zero_var=var_mat
  #var_mat_zero_exp_genes_set_to_zero_var[gene_w_exr_values==0]=NA
  melted_var_mat_zero_exp_genes_set_to_zero_var=melt(var_mat_zero_exp_genes_set_to_zero_var)
  colnames(melted_var_mat_zero_exp_genes_set_to_zero_var)[3]="Mean_Dev"
  #MEAN
  #mean_dev_by_taxon_exc_zero_exp_genes=cbind.data.frame(HMP32Name=as.character(var_mat_zero_exp_genes_set_to_zero_var$HMP32Name), Mean_Dev=as.numeric(rowMeans(var_mat_zero_exp_genes_set_to_zero_var[,2:ncol(var_mat_zero_exp_genes_set_to_zero_var)], na.rm=T)))
  #MEDIAN
  mean_dev_by_taxon_exc_zero_exp_genes=cbind.data.frame(HMP32Name=as.character(var_mat_zero_exp_genes_set_to_zero_var$HMP32Name), Mean_Dev=as.numeric(rowMedians(var_mat_zero_exp_genes_set_to_zero_var[,2:ncol(var_mat_zero_exp_genes_set_to_zero_var)], na.rm=T)))
  #GEOMETRIC MEAN
  #mean_dev_by_taxon_exc_zero_exp_genes=cbind.data.frame(HMP32Name=as.character(var_mat_zero_exp_genes_set_to_zero_var$HMP32Name), Mean_Dev=as.numeric(apply(var_mat_zero_exp_genes_set_to_zero_var[,2:ncol(var_mat_zero_exp_genes_set_to_zero_var)], 1, function(x) exp(mean(log(x), na.rm=T)))))
  
  #Seed Weight BLUPs
  #Load in phenotypes 
  #phenos=read.table("/home/kak268/work_data/maize_phenotypes/Ames_BLUPs/Merged_Ames_traits", header=T)
  phenos=read.table(file="/home/kak268/work_data/maize_phenotypes/282_panzea_phenotypes/traitMatrix_maize282NAM_v15-130212_Earweight_minus_Cobweight_melt_Means_and_BLUPS.txt", sep="\t", header=T)
  
  
  #Load in name converter so names canbe converted to Ames names
  name_decoder=read.table(header=T,"/media/kak268/A_2TB_Internal/RNA/Expressions_quants/samtools_idxstats_counts_bwamem_against_Zea_mays.AGPv3.29.cdna.all_w_Trinity_exc_NonAngiospermae_exc_over_85pct_B73_hits.fa/count_mat/renaming_scratch_folder/orig_1960names_and_alternate_names.txt")  
  name_decoder_L3Tip=name_decoder[name_decoder$TissueWODate=="L3Tip",]
  ###Make the decoder unique
  unique_name_decoder_L3Tip=name_decoder_L3Tip[!duplicated(name_decoder_L3Tip$HMP32Name),]
  
  
  #Merge dev by taxon with name decoder and phenos
  #inc zero exp gene dev
  #mean_dev_by_taxon_w_names=merge(mean_dev_by_taxon, unique_name_decoder_L3Tip, by="HMP32Name")
  #exc zero exp gene dev
  mean_dev_by_taxon_w_names=merge(melted_var_mat_zero_exp_genes_set_to_zero_var, unique_name_decoder_L3Tip, by="HMP32Name")
  
  mean_dev_by_taxon_w_names_unique_combs=mean_dev_by_taxon_w_names[!duplicated(mean_dev_by_taxon_w_names[c("HMP32Name", "AmesName", "Mean_Dev", "RNA_TaxaName", "RawPhenotypeNames")]),]
  
  #orig#mean_dev_by_taxon_w_names_and_phenos_unique_combs=merge(x=mean_dev_by_taxon_w_names_unique_combs, y=phenos, by.x="AmesName", by.y="TaxaName") 
  mean_dev_by_taxon_w_names_and_phenos_unique_combs=merge(x=mean_dev_by_taxon_w_names_unique_combs, y=phenos, by.x="RawPhenotypeNames", by.y="TaxaName") 
  
  
  #orig plot(mean_dev_by_taxon_w_names_and_phenos_unique_combs$Mean_Dev, mean_dev_by_taxon_w_names_and_phenos_unique_combs$PHT_BLUP, cex=1.5, col="deeppink4", pch=16, xlab=xlab, ylab=ylab)
  #pdf(paste(plot_out_dir,xlab,"vs_mean_seed_wt_BLUP.pdf", sep=""),width=11,height=11,paper='special') 
  plot(log10(mean_dev_by_taxon_w_names_and_phenos_unique_combs$Mean_Dev), mean_dev_by_taxon_w_names_and_phenos_unique_combs$SeedWeightBLUPs, col="darkorange3", pch=16, xlab=xlab, ylab="Seed Weight BLUPs")
  #orig line   line=lm(mean_dev_by_taxon_w_names_and_phenos_unique_combs$PHT_BLUP~mean_dev_by_taxon_w_names_and_phenos_unique_combs$Mean_Dev)
  line=lm(mean_dev_by_taxon_w_names_and_phenos_unique_combs$SeedWeightBLUPs~log10(mean_dev_by_taxon_w_names_and_phenos_unique_combs$Mean_Dev))
  abline(line)
  linesummary=summary(line)
  #legend("topright", bty="n", legend=as.expression(c(bquote(R^2 == .(round(linesummary$r.squared, digits=4))), paste("P = ", bquote(.(round(pf(linesummary$fstatistic[1], linesummary$fstatistic[2], linesummary$fstatistic[3], lower.tail = FALSE), digits=15))), sep=""))))# p == .(pval1) .(pval2)))
  legend("topright", bty="n", legend=as.expression(c(bquote(r == .(round(sqrt(linesummary$r.squared), digits=4))), paste("P = ", bquote(.(round(pf(linesummary$fstatistic[1], linesummary$fstatistic[2], linesummary$fstatistic[3], lower.tail = FALSE), digits=15))), sep=""))))# p == .(pval1) .(pval2)))
  intercept=round(linesummary$coefficients[1,1], 2)
  slope=round(linesummary$coefficients[2,1], 2)
  sign=NULL
  if (slope>=0){sign="+"}
  legend("top", bty="n", legend=paste("y = ", intercept, sign, slope,"x", sep=""))
  dev.off()
  
  #PHTBLUP
  #Load in phenotypes 
  rm(phenos)
  phenos=read.table("/home/kak268/work_data/maize_phenotypes/Ames_BLUPs/Merged_Ames_traits", header=T)
  rm(mean_dev_by_taxon_w_names_and_phenos_unique_combs)
  mean_dev_by_taxon_w_names_and_phenos_unique_combs=merge(x=mean_dev_by_taxon_w_names_unique_combs, y=phenos, by.x="AmesName", by.y="V1") 
  
  rm(linesummary)
  rm(line)
  #pdf(paste(plot_out_dir,xlab, "vs_mean_PHT_BLUP.pdf", sep=""),width=11,height=11,paper='special') 
  plot(log10(mean_dev_by_taxon_w_names_and_phenos_unique_combs$Mean_Dev), mean_dev_by_taxon_w_names_and_phenos_unique_combs$PHT_BLUP, col="cyan4", pch=16, xlab=xlab, ylab="Plant Height BLUP")
  line=lm(mean_dev_by_taxon_w_names_and_phenos_unique_combs$PHT_BLUP~log10(mean_dev_by_taxon_w_names_and_phenos_unique_combs$Mean_Dev))
  abline(line)
  linesummary=summary(line)
  #legend("topright", bty="n", legend=as.expression(c(bquote(R^2 == .(round(linesummary$r.squared, digits=4))), paste("P = ", bquote(.(round(pf(linesummary$fstatistic[1], linesummary$fstatistic[2], linesummary$fstatistic[3], lower.tail = FALSE), digits=15))), sep=""))))# p == .(pval1) .(pval2)))
  legend("topright", bty="n", legend=as.expression(c(bquote(r == .(round(sqrt(linesummary$r.squared), digits=4))), paste("P = ", bquote(.(round(pf(linesummary$fstatistic[1], linesummary$fstatistic[2], linesummary$fstatistic[3], lower.tail = FALSE), digits=15))), sep=""))))# p == .(pval1) .(pval2)))
  intercept=round(linesummary$coefficients[1,1], 2)
  slope=round(linesummary$coefficients[2,1], 2)
  sign=NULL
  if (slope>=0){sign="+"}
  legend("top", bty="n", legend=paste("y = ", intercept, sign, slope,"x", sep=""))
  dev.off()
  
  #GCA_FROM_YIELD_BLUPs
  #Load in phenotypes 
  rm(phenos)
  phenos=read.table("/home/kak268/work_data/maize_phenotypes/Yield_BLUPs_on_Testers/Ames_averaged_hyb_Yield_BLUPs_GCA.txt", header=T)
  phenos=phenos[complete.cases(phenos),]
  rm(mean_dev_by_taxon_w_names_and_phenos_unique_combs)
  mean_dev_by_taxon_w_names_and_phenos_unique_combs=merge(x=mean_dev_by_taxon_w_names_unique_combs, y=phenos, by.x="AmesName", by.y="AmesName") 
  mean_dev_by_taxon_w_names_and_phenos_unique_combs_by_RNAname=merge(x=mean_dev_by_taxon_w_names_unique_combs, y=phenos, by.x="RNA_TaxaName" , by.y="AmesName")
  mean_dev_by_taxon_w_names_and_phenos_unique_combs_total=rbind(mean_dev_by_taxon_w_names_and_phenos_unique_combs_by_RNAname, mean_dev_by_taxon_w_names_and_phenos_unique_combs)
  print(dim(mean_dev_by_taxon_w_names_and_phenos_unique_combs))
  print(dim(mean_dev_by_taxon_w_names_and_phenos_unique_combs_by_RNAname))
  print(dim(mean_dev_by_taxon_w_names_and_phenos_unique_combs_total))
  #print(mean_dev_by_taxon_w_names_and_phenos_unique_combs_total)
  
  rm(linesummary)
  rm(line)
  #pdf(paste(plot_out_dir,xlab,"vs_mean_yield_BLUP.pdf", sep=""),width=11,height=11,paper='special') 
  plot(log10(mean_dev_by_taxon_w_names_and_phenos_unique_combs_total$Mean_Dev), mean_dev_by_taxon_w_names_and_phenos_unique_combs_total$Mean_yield_BLUP, col="darkred", pch=16, xlab=xlab, ylab="Mean Yield BLUP")
  line=lm(mean_dev_by_taxon_w_names_and_phenos_unique_combs_total$Mean_yield_BLUP~log10(mean_dev_by_taxon_w_names_and_phenos_unique_combs_total$Mean_Dev))
  abline(line)
  linesummary=summary(line)
  #legend("topright", bty="n", legend=as.expression(c(bquote(R^2 == .(round(linesummary$r.squared, digits=4))), paste("P = ", bquote(.(round(pf(linesummary$fstatistic[1], linesummary$fstatistic[2], linesummary$fstatistic[3], lower.tail = FALSE), digits=15))), sep=""))))# p == .(pval1) .(pval2)))
  legend("topright", bty="n", legend=as.expression(c(bquote(r == .(round(sqrt(linesummary$r.squared), digits=4))), paste("P = ", bquote(.(round(pf(linesummary$fstatistic[1], linesummary$fstatistic[2], linesummary$fstatistic[3], lower.tail = FALSE), digits=15))), sep=""))))# p == .(pval1) .(pval2)))
  intercept=round(linesummary$coefficients[1,1], 2)
  slope=round(linesummary$coefficients[2,1], 2)
  sign=NULL
  if (slope>=0){sign="+"}
  legend("top", bty="n", legend=paste("y = ", intercept, sign, slope,"x", sep=""))
  dev.off()
}




phenos=read.table("/home/kak268/work_data/maize_phenotypes/Yield_BLUPs_on_Testers/Ames_averaged_hyb_Yield_BLUPs_GCA.txt", header=T)




#Random 5000 genes
basedir="/media/kak268/A_2TB_Internal/RNA/Expressions_quants/HTSEQ_counts_STAR_against_Zea_mays_B73_AGPv3.29_w_gtf_annotation/count_mat/count_mats_matched_to_AltNames/"
rand_5000_gene_w_exr_values=read.table(paste(basedir, "L3tip_vs_rare_alleles/L3Tip_rand_5000_expressed_genes_w_STAR_HTSEQ_DESEQ_actual_counts.txt", sep=""), header=T)
#rand_5000_gene_w_exr_values=read.table(paste(basedir, "L3tip_vs_rare_alleles/L3Tip_rand_5000_expressed_genes_w_STAR_HTSEQ_DESEQ_actual_counts.txt", sep=""), header=T)

plot_exp_dev_v_phenos(rand_5000_gene_w_exr_values, xlab="Mean sq. deviation of individual expression from population mean in 5000 rand genes")
plot_exp_dev_v_phenos(rand_5000_gene_w_exr_values, xlab="log10 Mean sq. deviation of individual expression from population mean in 5000 rand genes")





#TEST of plotting and reshaping dataframes

df=cbind(a=c("a","b","c","d", "d", "d", "e", "e"), x=1:8, y=3:11)
df=cbind(v=rnorm(n=100,mean=0,sd=5), w=rnorm(n=100,mean=10,sd=5),x=rnorm(n=100,mean=20,sd=5),y=rnorm(n=100,mean=30,sd=5),z=rnorm(n=100,mean=40,sd=5)
df2=cbind.data.frame(a=c(rep(1, 10), rep(2, 10), rep(3, 10), rep(4, 10), rep(5, 10)), b=c(rnorm(n=10,mean=10,sd=5),rnorm(n=10,mean=20,sd=5),rnorm(n=10,mean=30,sd=5),rnorm(n=10,mean=40,sd=5), rnorm(n=10,mean=50,sd=5)), c=c(rep(1.1, 10), rep(2, 10), rep(3, 10), rep(4, 10), rep(5, 10)))
df2[1,2]=NA


#library(reshape)
#meltd_df=melt(df)
#plot(meltd_df$X2, meltd_df$value)
plot(df2$a, df2$b)
line=lm(df2$b~df2$a)
summary(line)
summary(line)$p-value
summary(line)$r.squared
abline(line)

df3=data.frame(a=double(), b=double())
for (i in unique(df2$a)){
  df3[i,1]=i
  df3[i,2]=mean(df2[df2$a==i,2])
}
plot(df3$a, df3$b)
line=lm(df3$b~df3$a)
summary(line)
summary(line)$r.squared
abline(line)
