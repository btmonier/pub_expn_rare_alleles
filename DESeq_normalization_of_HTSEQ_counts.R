# Karl Kremling
#Jan 23 2016
#DESeq norm of samtools idxstats counts

#source("http://bioconductor.org/biocLite.R")
#import DESeq and normalize
#biocLite("DESeq2")
# if this doesn't work on linux you amy need to sudo apt-get install libxml2-dev
library("DESeq2")

#install.packages("fastmatch")
library("fastmatch")

library("parallel")

#AGP v3.29 version 
basedir="/media/kak268/A_2TB_Internal/RNA/Expressions_quants/HTSEQ_counts_STAR_against_Zea_mays_B73_AGPv3.29_w_gtf_annotation/"
#AGP v4.34 version w unique mappers only (but shouldn't matter because HTSeq only uses top one)
#basedir="/media/kak268/A_2TB_Internal/RNA/Expressions_quants/HTSEQ_counts_STAR_against_Zea_mays_B73_AGPv4.34_w_gtf_annotation/HTSEQ_counts_multimap1/"
#AGP v4.34 version w up to 10 multimapping (but shouldn't matter because HTSeq only uses top one)
#basedir="/media/kak268/A_2TB_Internal/RNA/Expressions_quants/HTSEQ_counts_STAR_against_Zea_mays_B73_AGPv4.34_w_gtf_annotation/HTSEQ_counts_multimap10/"

#AGP v4.34 version 3' end extended 500 bp
#basedir="/media/kak268/A_2TB_Internal/RNA/Expressions_quants/HTSEQ_counts_STAR_against_Zea_mays_B73_AGPv4.34_w_gtf_annotation_genes_ext_500bp/"

#Sorghum Btx623 base and tip
#basedir="/media/kak268/A_2TB_Internal/RNA/Expressions_quants/HTSEQ_counts_STAR_against_sorghum_and_maize/all_tissues_08262017sorghum_3p_FWD/HTSEQ_counts/sorghum_only/"


addNoise <- function(mtx) {
  set.seed(122)
  set.seed(123) # have to do this to reset it to the beignning , also wasnt using set seed for the individual tissues
  a=mtx[(!is.na(mtx))]
  max_val=min(a[a>0])/2
  temp_mtx = mtx
  cat(max_val)# determines the min non-zero expression value which will be the upper bound for the random number generation
  if (!is.matrix(mtx)) mtx <- matrix(mtx, byrow = TRUE, nrow = 1) # converts the df to a matrix
  random.stuff <- matrix(runif(prod(dim(mtx)), min = 0.00000001, max = max_val), nrow = dim(mtx)[1]) #creates a matrix of random values 
  random.stuff + temp_mtx
}

# addNoise <- function(mtx) {
#   if (!is.matrix(mtx)) mtx <- matrix(mtx, byrow = TRUE, nrow = 1)
#   random.stuff <- matrix(runif(prod(dim(mtx)), min = 0.00000001, max = 0.0001), nrow = dim(mtx)[1])
#   cat(length(random.stuff), "\n")
#   cat(length(mtx))
#   random.stuff + mtx
# }


###
#Function to take in HTSEQ counts and convert them into a normalized matrix
#cts_to_norm_mat <- function(basedir, tissue){
date()
files <- list.files(paste(basedir, "HTSEQ_counts/", sep = ""), full.names=T, pattern="count.txt", recursive=F)
#get rid of the empty files
non0list=NULL
for (i in files){non0list=c(non0list, file.info(i)$size!=0)}
files=files[non0list]
genes <- read.table(files[1], header=F, sep="\t")[,1]    # gene/contig names MAKE SURE ALL FILES HAVE SAME GENE ORDER
df    <- do.call(cbind,mclapply(files,function(fn)read.table(fn,header=F, sep="\t")[,2]))  # takes 30 mins for 2208 files w 150k genes
#colnames(df)= substring(files, first=)
trim_left=sub(".*BGXX_", "", files)
trim_right=sub("_R1_.*", "", trim_left)
colnames(df) = trim_right
rownames(df)= genes  
date()

#Remove the complete count columns at the end of the HTSEQ vectors
df=df[!rownames(df) %in% c("__no_feature","__ambiguous", "__too_low_aQual", "__not_aligned", "__alignment_not_unique"),]


###Use the poisitive id list from the genetic distances to keep only the counts which are from positively identified samples 
positive_id_file_list=read.table("/media/kak268/B_2TB_Internal/Genotypes/RNA_SNPs_Chr10_Fei_pipeline/Positive_id_list_filenames_only_RNAseq_samples_based_on_dist_from_corresponding_HMP32_line.txt")
trim_left_positive_list=sub(".*BGXX_", "", positive_id_file_list[,1])
trim_right_positive_list=sub("_R1.fastq.gz", "", trim_left_positive_list)
#keep only cols with positive id
df_matched_by_genet_dist=subset(df, select=trim_right_positive_list)
#dim(df_redundant_trinity_summed_and_orig_B73_counts_only_positively_matched_from_genet_dist)
#AGP v3.29 version
#write.table(df_matched_by_genet_dist, "/media/kak268/A_2TB_Internal/RNA/Expressions_quants/HTSEQ_counts_STAR_against_Zea_mays_B73_AGPv3.29_w_gtf_annotation/count_mat/df_STAR_HTSeq_counts_B73.txt", quote=F) 
#AGP v4.34 version multimap1
#write.table(df_matched_by_genet_dist, "/media/kak268/A_2TB_Internal/RNA/Expressions_quants/HTSEQ_counts_STAR_against_Zea_mays_B73_AGPv4.34_w_gtf_annotation/count_mat/df_STAR_HTSeq_counts_B73.txt", quote=F) 
#AGP v4.34 version multimap10
#write.table(df_matched_by_genet_dist, "/media/kak268/A_2TB_Internal/RNA/Expressions_quants/HTSEQ_counts_STAR_against_Zea_mays_B73_AGPv4.34_w_gtf_annotation/count_mat_multimap10/df_STAR_HTSeq_counts_B73.txt", quote=F) 

##AGP v4.34 version 3' end extended 500 bp
#write.table(df_matched_by_genet_dist, paste(basedir, "count_mat/df_STAR_HTSeq_counts_B73.txt", sep = ""), quote=F) 





###DESEQ Normalization

#Make design matrix for use in DESeq2
colData=colnames(df_matched_by_genet_dist)
#colnames(colData) = c("Individual")


#Make improved taxanames
#TaxaNames=sub(paste(".*", tissue, "_", sep=""), "", colData[,1]) #remove everything up to and including the tissue name ".
#TaxaNames=substr(TaxaNames, 1, nchar(TaxaNames)-7)

#Deseq normalization
dds =  DESeqDataSetFromMatrix(countData  = df_matched_by_genet_dist,  colData  = as.data.frame(colData), design= ~ 1) # this colData and design are irrelevant since all we want are counts normed by total numbers, not by design
dds = estimateSizeFactors(dds)
counts.mat = counts(dds, normalized=T)

#Remove rows (genes/contigs) which have 0 expression values in all the taxa (or > 3/4 of the taxa)
#counts.mat = counts.mat[apply(counts.mat==0,1,sum)<=0.75*dim(counts.mat)[2],]
counts.mat = counts.mat[apply(counts.mat==0,1,sum)<=dim(counts.mat)[2]-1,]
counts.mat = counts.mat[ order(row.names(counts.mat)), ] # sort by the name of the gene so that specific genes can be easily compared between df
#write counts without rounding, or adding small rand and log2 transforming
counts.mat.w.smp.names=rbind(colnames(df_matched_by_genet_dist), counts.mat)
rownames(counts.mat.w.smp.names)[1]="<Trait>" # this allows it to be recognized by TASSEL. This works because write.table puts first col name over col 1, which results in an off by 1 error unless you add a string there, or include the option col.names=NA when you write to file
#AGP v3.29 version
write.table(x=t(counts.mat.w.smp.names), file="/media/kak268/A_2TB_Internal/RNA/Expressions_quants/HTSEQ_counts_STAR_against_Zea_mays_B73_AGPv3.29_w_gtf_annotation/count_mat/df_STAR_HTSeq_counts_B73_match_based_on_genet_dist_DESeq2_normed.txt", quote=F, col.names=T, row.names=F)
#AGP v4.34 version
#write.table(x=t(counts.mat.w.smp.names), file=paste(basedir, "/count_mat/df_STAR_HTSeq_counts_B73_match_based_on_genet_dist_DESeq2_normed.txt",sep=""), quote=F, col.names=T, row.names=F)

#DESeq normalization and FPM conversion
fpm.counts.mat=fpm(dds, robust = TRUE)
fpm.counts.mat = fpm.counts.mat[apply(fpm.counts.mat==0,1,sum)<=dim(fpm.counts.mat)[2]-1,]
fpm.counts.mat = fpm.counts.mat[ order(row.names(fpm.counts.mat)), ] # sort by the name of the gene so that specific genes can be easily compared between df
#write counts without rounding, or adding small rand and log2 transforming
fpm.counts.mat.w.smp.names=rbind(colnames(df_matched_by_genet_dist), fpm.counts.mat)
rownames(fpm.counts.mat.w.smp.names)[1]="<Trait>" # this allows it to be recognized by TASSEL. This works because write.table puts first col name over col 1, which results in an off by 1 error unless you add a string there, or include the option col.names=NA when you write to file
#AGP v3.29 version
write.table(x=t(fpm.counts.mat.w.smp.names), file="/media/kak268/A_2TB_Internal/RNA/Expressions_quants/HTSEQ_counts_STAR_against_Zea_mays_B73_AGPv3.29_w_gtf_annotation/count_mat/df_STAR_HTSeq_counts_B73_match_based_on_genet_dist_DESeq2_normed_fpm.txt", quote=F, col.names=T, row.names=F)
#AGP v4.34 version
#write.table(x=t(fpm.counts.mat.w.smp.names), file=paste(basedir, "/count_mat/df_STAR_HTSeq_counts_B73_match_based_on_genet_dist_DESeq2_normed_fpm.txt",sep=""), quote=F, col.names=T, row.names=F)



#write counts without rounding, but add small rand and log2 transform
counts.mat.w.sm.add.log2=log2(addNoise(counts.mat))
counts.mat.w.sm.add.log2.w.smp.names=rbind(colnames(df_matched_by_genet_dist), counts.mat.w.sm.add.log2)
rownames(counts.mat.w.sm.add.log2.w.smp.names)[1]="<Trait>" # this allows it to be recognized by TASSEL. This works because write.table puts first col name over col 1, which results in an off by 1 error unless you add a string there, or include the option col.names=NA when you write to file
#AGP v3.29 version
#write.table(x=t(counts.mat.w.sm.add.log2.w.smp.names), file="/media/kak268/A_2TB_Internal/RNA/Expressions_quants/HTSEQ_counts_STAR_against_Zea_mays_B73_AGPv3.29_w_gtf_annotation/count_mat/df_STAR_HTSeq_counts_B73_match_based_on_genet_dist_DESeq2_normed_sm_rand_add_log2tform.txt", quote=F, col.names=T, row.names=F)
#AGP v4.34 version
#write.table(x=t(counts.mat.w.sm.add.log2.w.smp.names), file=paste(basedir, "/count_mat/df_STAR_HTSeq_counts_B73_match_based_on_genet_dist_DESeq2_normed_sm_rand_add_log2tform.txt", sep=""), quote=F, col.names=T, row.names=F)


#write DESeq2 counts with rounding, but without adding small rand and log2 transforming
counts.mat.round=round(counts.mat,digits=5)
counts.mat.round.w.smp.names=rbind(colnames(df_matched_by_genet_dist), counts.mat.round)
rownames(counts.mat.round.w.smp.names)[1]="<Trait>" # this allows it to be recognized by TASSEL. This works because write.table puts first col name over col 1, which results in an off by 1 error unless you add a string there, or include the option col.names=NA when you write to file
#AGP v3.29 version
write.table(x=t(counts.mat.round.w.smp.names), file="/media/kak268/A_2TB_Internal/RNA/Expressions_quants/HTSEQ_counts_STAR_against_Zea_mays_B73_AGPv3.29_w_gtf_annotation/count_mat/df_STAR_HTSeq_counts_B73_match_based_on_genet_dist_DESeq2_normed_rounded.txt", quote=F, col.names=T, row.names=F)
#AGP v4.34 version  w unique mappers
#write.table(x=t(counts.mat.round.w.smp.names), file=paste(basedir, "/count_mat/df_STAR_HTSeq_counts_B73_match_based_on_genet_dist_DESeq2_normed_rounded.txt",sep=""), quote=F, col.names=T, row.names=F)

#write FPM DESeq2 counts with rounding, but without adding small rand and log2 transforming
fpm.counts.mat.round=round(fpm.counts.mat,digits=5)
fpm.counts.mat.round.w.smp.names=rbind(colnames(df_matched_by_genet_dist), fpm.counts.mat.round)
rownames(fpm.counts.mat.round.w.smp.names)[1]="<Trait>" # this allows it to be recognized by TASSEL. This works because write.table puts first col name over col 1, which results in an off by 1 error unless you add a string there, or include the option col.names=NA when you write to file
#AGP v3.29 version
write.table(x=t(fpm.counts.mat.round.w.smp.names), file="/media/kak268/A_2TB_Internal/RNA/Expressions_quants/HTSEQ_counts_STAR_against_Zea_mays_B73_AGPv3.29_w_gtf_annotation/count_mat/df_STAR_HTSeq_counts_B73_match_based_on_genet_dist_DESeq2_normed_fpm_rounded.txt", quote=F, col.names=T, row.names=F)
#AGP v4.34 version  w unique mappers
#write.table(x=t(counts.mat.round.w.smp.names), file=paste(basedir, "/count_mat/df_STAR_HTSeq_counts_B73_match_based_on_genet_dist_DESeq2_normed_rounded.txt",sep=""), quote=F, col.names=T, row.names=F)



#write counts with rounding after adding small rand and log2 transform
counts.mat.w.sm.add.log2.rd=round(log2(addNoise(counts.mat)), digits=5)
counts.mat.w.sm.add.log2.rd.w.smp.names=rbind(colnames(df_matched_by_genet_dist), counts.mat.w.sm.add.log2.rd)
rownames(counts.mat.w.sm.add.log2.rd.w.smp.names)[1]="<Trait>" # this allows it to be recognized by TASSEL. This works because write.table puts first col name over col 1, which results in an off by 1 error unless you add a string there, or include the option col.names=NA when you write to file
#AGP v3.29 version
#write.table(x=t(counts.mat.w.sm.add.log2.rd.w.smp.names), file="/media/kak268/A_2TB_Internal/RNA/Expressions_quants/HTSEQ_counts_STAR_against_Zea_mays_B73_AGPv3.29_w_gtf_annotation/count_mat/df_STAR_HTSeq_counts_B73_match_based_on_genet_dist_DESeq2_normed_sm_rand_add_log2tform_rounded.txt", quote=F, col.names=T, row.names=F)
#AGP v4.34 version
write.table(x=t(counts.mat.w.sm.add.log2.rd.w.smp.names), file=paste(basedir, "/count_mat/df_STAR_HTSeq_counts_B73_match_based_on_genet_dist_DESeq2_normed_sm_rand_add_log2tform_rounded.txt",sep=""), quote=F, col.names=T, row.names=F)



#read in DeSEQ counts write counts by tissue after removing genes that are completely zero in a tissue adding small rand and boxcox transforming
filename="df_STAR_HTSeq_counts_B73_match_based_on_genet_dist_DESeq2_normed_rounded_origNames_and_Altnames_noL3Mid" #exclude .txt
#AGP v3.29 version
count_dir="/media/kak268/A_2TB_Internal/RNA/Expressions_quants/HTSEQ_counts_STAR_against_Zea_mays_B73_AGPv3.29_w_gtf_annotation/count_mat/count_mats_matched_to_AltNames/"
#AGP v4.34 version
count_dir="/media/kak268/A_2TB_Internal/RNA/Expressions_quants/HTSEQ_counts_STAR_against_Zea_mays_B73_AGPv4.34_w_gtf_annotation/count_mat/count_mats_matched_to_AltNames/"

DESeq_counts=read.table(file=paste(count_dir, filename, ".txt", sep=""), header=T) # takes ~5 mins for 2000x140000
DESeq_counts[DESeq_counts=="NA"]=NaN
DESeq_counts=DESeq_counts
last_name_col=10
library("MASS")
for (tissue in unique(DESeq_counts$TissueWODate)){
  cat(date(), " ")
  cat(tissue, "\n")
  subset_DESeq_counts=DESeq_counts[DESeq_counts$TissueWODate==tissue,] #subset by tissue
  subset_names=subset_DESeq_counts[,1:last_name_col]
  subset_no_names=subset_DESeq_counts[,-c(1:last_name_col)]
  subset_no_names_no_complete_0 = subset_no_names[,apply(subset_no_names==0,2,sum)<=dim(subset_no_names)[1]-1] # remove cols with only zeroes (for cols make sure to use 2 in apply and put the comma before the expression)
  subset_no_names_sm_rand=addNoise(subset_no_names_no_complete_0) #adds a small random value which is < 1/2 the minimum expression value in the df
  lambdaMASS_vec=NULL
  for (i in colnames(subset_no_names_sm_rand)){
    bc= boxcox(subset_no_names_sm_rand[[i]]~1, plotit=F, lambda=seq(-2,2,0.01))
    lambdaMASS_vec<-c(lambdaMASS_vec, bc$x[which.max(bc$y)])
  }
  #http://www.isixsigma.com/tools-templates/normality/making-data-normal-using-box-cox-power-transformation/
  bc_tranformer = function(obs, lambda) {
    y=NULL
    if (lambda!=0){y=obs^lambda}
    else {y=log(obs)}
    return(y)
  }
  df_bc=mapply(bc_tranformer, subset_no_names_sm_rand, lambdaMASS_vec)
  df_bc_w_names=cbind(subset_names, df_bc)
  write.table(x=df_bc_w_names, file=paste(count_dir, "count_mats_separated_by_tissue_and_individually_box_coxed/", filename, "_sm_rand_and_box_coxed_", tissue,".txt", sep=""), row.names=F, quote=F)
}
test=cbind(x=c(1,2,3,4,5,0), y=c(3,4,11,6,0.1,7))
addNoise(test)



#Apply test and box cox test with lambda=0
# subset_no_names_sm_rand=subset_no_names_sm_rand[1:10, 1:10]
# lambdaMASS_vec=NULL
# for (i in colnames(subset_no_names_sm_rand)){
#   bc= boxcox(subset_no_names_sm_rand[[i]]~1, plotit=F, lambda=seq(-2,2,0.01))
#   lambdaMASS_vec<-c(lambdaMASS_vec, bc$x[which.max(bc$y)])
# }
# lambdaMASS_vec[5]=0 # hardcoding one of the lamdas to 0 in order to test the tranformation code
# bc_tranformer = function(obs, lambda) {
#   y=NULL
#   if (lambda!=0){y=obs^lambda}
#   else {y=log(obs)}
#   return(y)
# }
# df_bc=mapply(bc_tranformer, subset_no_names_sm_rand, lambdaMASS_vec)
# df_bc[1:2, 1:5]
