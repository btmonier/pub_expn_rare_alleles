#Karl Kremling
#filter top eQTL hits from TASSEl into CIS only
#install.packages('readr')
library(readr)
#install.packages('dplyr')
library(dplyr)
#install.packages("plyr")
library(plyr)
#install.packages("data.table")
library(data.table)

cis_radius=1000000
eqtl_dir="/media/kak268/A_2TB_Internal/RNA/eQTL_associations/eQTL_associations_avgDups/results/"

#pull genes from gff file and keep only the genes which are on chr1 - chr10
gff_file=read.delim("/home/kak268/work_data/maize_genome/Zea_mays_gene_only.AGPv3.29.gff3", header=F, sep="\t", quote="", stringsAsFactors = F) # add skip=6 if it's a full gff file
#gff_file=read.delim("/local/workdir/kak268/rare_allel_calc/Zea_mays_gene_only.AGPv3.29.gff3", header=F, sep="\t", quote="") # add skip=6 if it's a full gff file

colnames(gff_file)=c("Chr", "Annot_source", "Feature_type", "Start", "End", "blankcol1", "Strand", "blankcol2", "Feature_Name")
transcript_file=gff_file[gff_file$Feature_type=='gene',] # keep only the entries in the GFF file which are transcripts/genes
transcript_file$Feature_Name=sub(pattern= "ID=gene:", replacement='', x=transcript_file$Feature_Name)
transcript_file$Feature_Name=sub(pattern= ";.*", replacement='', x=transcript_file$Feature_Name)
transcript_file_only_chr1_to_chr10=transcript_file[transcript_file$Chr %in% c(1,2,3,4,5,6,7,8,9,10),] # keeps all transcript file info
transcript_file_only_chr1_to_chr10_gene_list=transcript_file[transcript_file$Chr %in% c(1,2,3,4,5,6,7,8,9,10),"Feature_Name"] # keeps only the genes on one of the ten chromosomes

colnames(transcript_file_only_chr1_to_chr10)[4:5]=c("Gene.Start", "Gene.End")
colnames(transcript_file_only_chr1_to_chr10)[1]=c("Gene.Chr")
transcript_file_only_chr1_to_chr10$Gene.Chr=as.numeric(transcript_file_only_chr1_to_chr10$Gene.Chr)

colnames(transcript_file_only_chr1_to_chr10)[9]="Trait"

# temp=tail(transcript_file_only_chr1_to_chr10)
# temp2=temp
# temp2$Gene.Chr=as.numeric(temp$Gene.Chr)

for (tissue in c("L3Tip")){#, "L3Base", "Kern", "LMAD", "LMAN", "GRoot", "GShoot")){
  all_chrs=NULL
  chr=NULL
  for (chr in 1:10){
    cat("tissue and chr are ", tissue, " ", chr)
    eqtl_results=read_delim(paste(eqtl_dir,"STAR_HTSEQ_DESEQ_BoxCox_sm_rand_25_PEERs_5MDSPC_from_all/",tissue,"_merged_flt_c", chr,"_maxP000001.hmp321.onlyRNAset_MAFover005.KNNi.w.df_STAR_HTSeq_counts_B73_match_based_on_genet_dist_DESeq2_normed_rounded_origNames_and_Altnames_sm_rand_and_box_coxed_w_MDSPCs_and_25PEERs_avgDups.txt", sep = ""), delim = "\t", col_names = T)

    eqtl_results_w_gff=inner_join(eqtl_results, transcript_file_only_chr1_to_chr10, by="Trait")
    eqtl_results_w_gff_sing_chr=eqtl_results_w_gff[eqtl_results_w_gff$Gene.Chr==chr,]
    #eqtl_results_w_gff_sing_chr_temp=eqtl_results_w_gff_sing_chr
    #eqtl_results_w_gff_sing_chr=eqtl_results_w_gff_sing_chr[1:10,]                         
    eqtl_results_w_gff_sing_chr_cis=eqtl_results_w_gff_sing_chr[abs(eqtl_results_w_gff_sing_chr$Pos-eqtl_results_w_gff_sing_chr$Gene.Start)<=cis_radius | abs(eqtl_results_w_gff_sing_chr$Pos-eqtl_results_w_gff_sing_chr$Gene.End)<=cis_radius,]
    
    #sort by SNP name and secondarily r2 within SNP name
    eqtl_results_w_gff_sing_chr_cis_sorted=eqtl_results_w_gff_sing_chr_cis[order(eqtl_results_w_gff_sing_chr_cis$Pos, eqtl_results_w_gff_sing_chr_cis$r2, decreasing = T),]
    #keep only the first entry per SNP, which is the one with the highest R2
    eqtl_results_w_gff_sing_chr_cis_sorted_unique=eqtl_results_w_gff_sing_chr_cis_sorted[!duplicated(eqtl_results_w_gff_sing_chr_cis_sorted$Marker),]
    #eqtl_results_w_gff_sing_chr_cis_unique_SNPs=ddply(eqtl_results_w_gff_sing_chr_cis, "Pos", function(z) max(z$r2)) 
    eqtl_results_w_gff_sing_chr_cis_sorted_unique=eqtl_results_w_gff_sing_chr_cis_sorted_unique[,1:7]
    if(!exists("all_chrs")){
      all_chrs=eqtl_results_w_gff_sing_chr_cis_sorted_unique
    }else{
      all_chrs=rbind(all_chrs,eqtl_results_w_gff_sing_chr_cis_sorted_unique)
    }
  }
  write.table(all_chrs, file = paste(eqtl_dir, "top_cis_hit_dir/","top_hit_by_SNP_",tissue,"_merged_flt_allchr_maxP000001.hmp321.onlyRNAset_MAFover005.KNNi.w.df_STAR_HTSeq_counts_B73_match_based_on_genet_dist_DESeq2_normed_rounded_origNames_and_Altnames_sm_rand_and_box_coxed_w_MDSPCs_and_25PEERs_avgDups.txt", sep = ""), quote=F, sep = "\t", row.names = F)
}

tissue="LMAN"

a=c(1,2,3,6,78)
b=c(2,1,5,4,90)

c=cbind(a,b)
c.dt=as.data.table(c)
pmin(c.dt)
