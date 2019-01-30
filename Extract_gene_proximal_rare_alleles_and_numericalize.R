#Karl Kremling
#Aril 20 2016
#Make list of SNPs near genes, extract these sites from hmp file of rare alleles, make numerical genotype file of these gene-proximal rare alleles 

#install.packages("miscTools")
library("miscTools")
library("readr")
#install.packages('reshape')
library("reshape")
#install.packages('doParallel')
library("doParallel")
#install.packages("registerDoMc")
#library("registerDoMc")

#install.packages("doMC")
library("doMC")
registerDoMC(cores=12)
#
#
#
###Read in gene exp levels
#NOTE THIS IS FROM UN TFORMED GENES, only DESEQ-ed
basedir="/media/kak268/A_2TB_Internal/RNA/Expressions_quants/HTSEQ_counts_STAR_against_Zea_mays_B73_AGPv3.29_w_gtf_annotation/count_mat/count_mats_matched_to_AltNames/"
expression_df=read.table(paste(basedir, "df_STAR_HTSeq_counts_B73_match_based_on_genet_dist_DESeq2_normed_rounded_origNames_and_Altnames.txt", sep=""), header=T)


#pull genes from gff file and keep only the genes which are on chr1 - chr10
gff_file=read.delim("/home/kak268/work_data/maize_genome/Zea_mays_gene_only.AGPv3.29.gff3", header=F, sep="\t", quote="") # add skip=6 if it's a full gff file
#gff_file=read.delim("/local/workdir/kak268/rare_allel_calc/Zea_mays_gene_only.AGPv3.29.gff3", header=F, sep="\t", quote="") # add skip=6 if it's a full gff file

colnames(gff_file)=c("Chr", "Annot_source", "Feature_type", "Start", "End", "blankcol1", "Strand", "blankcol2", "Feature_Name")
transcript_file=gff_file[gff_file$Feature_type=='gene',] # keep only the entries in the GFF file which are transcripts/genes
transcript_file$Feature_Name=sub(pattern= "ID=gene:", replacement='', x=transcript_file$Feature_Name)
transcript_file$Feature_Name=sub(pattern= ";.*", replacement='', x=transcript_file$Feature_Name)
transcript_file_only_chr1_to_chr10=transcript_file[transcript_file$Chr %in% c(1,2,3,4,5,6,7,8,9,10),] # keeps all transcript file info
transcript_file_only_chr1_to_chr10_gene_list=transcript_file[transcript_file$Chr %in% c(1,2,3,4,5,6,7,8,9,10),"Feature_Name"] # keeps only the genes on one of the ten chromosomes



#read in HMP file to be used by each tissue
#All SNPs up to 005 MAF for all RNAset
HMP_file=read_delim("/media/kak268/B_2TB_Internal/Genotypes/HM32/LDKNNi/SNPNames_merged_flt_c1-10.hmp321.onlyRNAset_MAFover0under005.KNNi.hmp.txt", delim=" ", col_names=T)


#make individual chrs
for (i in unique(HMP_file$chrom)){
  cat(date())
  assign(paste("HMP_file_chr", i, sep = ""), HMP_file[HMP_file$chrom==i,])
}
rm(HMP_file)
rm(i)
gc(reset=T)

SNPid_list=NULL
features=unique(transcript_file_only_chr1_to_chr10$Feature_Name)
#SNPid_list1=foreach(i=1:20000, .combine=rbind) %dopar% {
SNPid_list=foreach(i=1:(length(features)), .combine=rbind) %dopar% {
  
  line=transcript_file_only_chr1_to_chr10[transcript_file_only_chr1_to_chr10$Feature_Name==features[i],]
  cat(i, "\n")
  print(line)
  #SNPs=NULL
  #cat(date())
  #HMP_file_single_chr=HMP_file[HMP_file$chrom==line$Chr,]
  #cat(date())
  
  HMP_file_single_chr=eval(parse(text=(paste("HMP_file_chr", line$Chr, sep = ""))))
  print(head(HMP_file_single_chr))
  if (line$Strand=="-"){
    cat("in - if \n")
    fivePend=line$End
    data.frame(HMP_file_single_chr[c(HMP_file_single_chr$pos<=(fivePend+5000) & HMP_file_single_chr$pos>fivePend),]) # getting the SNPs 5kb upstream
  } else if (line$Strand=="+"){
    cat("in + if \n")
    fivePend=line$Start
    data.frame(HMP_file_single_chr[c(HMP_file_single_chr$pos>=(fivePend-5000) & HMP_file_single_chr$pos<fivePend),])
  }
}
Sys.time()

#SNPid_list=rbind(SNPid_list1, SNPid_list2)

date()
snp_dir="/media/kak268/A_2TB_Internal/RNA/Expressions_quants/HTSEQ_counts_STAR_against_Zea_mays_B73_AGPv3.29_w_gtf_annotation/count_mat/count_mats_matched_to_AltNames/expr_vs_rare_alleles_AvgOverDups/Upstream_rare_allele_counts_with_expr_bins/all_genes/"
#write.table(SNPid_list$rs., paste(snp_dir, "5kb_upstream_SNPlist_",h,"_MAFover0under005_top_5000_expressed_genes_w_STAR_HTSEQ_DESEQ_median_counts.txt", sep=""), quote=F, row.names=F)
write.table(SNPid_list$rs., paste(snp_dir, "5kb_upstream_SNPlist_MAFover0under005_all_genes_w_STAR_HTSEQ_DESEQ.txt", sep=""), quote=F, row.names=F)

#Last checked with TASSEL 5.2.40 also works and gives identical md5
#Extract sites based on w/in 5kb of gene site list made above
string1=paste("perl /home/kak268/Software/tassel-5-standalone/run_pipeline.pl -Xmx50g -fork1 -h /media/kak268/B_2TB_Internal/Genotypes/HM32/LDKNNi/merged_flt_c1-10.hmp321.onlyRNAset_MAFover0under005.KNNi.hmp.txt -filterAlign -includeSiteNamesInFile ",snp_dir, "5kb_upstream_SNPlist_MAFover0under005_all_genes_w_STAR_HTSEQ_DESEQ.txt -export ",snp_dir, "merged_flt_c1-10.hmp321.RNAsetMAFover0under005_w_in_5kb_of_allgenes.KNNi -runfork1", sep="")

#numericalize
string3=paste("perl /home/kak268/Software/tassel-5-standalone/run_pipeline.pl -Xmx50g -fork1 -h ", snp_dir, "merged_flt_c1-10.hmp321.RNAsetMAFover0under005_w_in_5kb_of_allgenes.KNNi.hmp.txt -NumericalGenotypePlugin -endPlugin -export ", snp_dir, "merged_flt_c1-10.hmp321.RNAsetMAFover0under005_w_in_5kb_of_allgenes.KNNi.hmp.numeric -exportType ReferenceProbability", sep="")

#transpose and convert 0.5 hets to 3 for reading by R
string5=paste("tail -n +2 ", snp_dir, "merged_flt_c1-10.hmp321.RNAsetMAFover0under005_w_in_5kb_of_allgenes.KNNi.hmp.numeric.txt | bash /media/kak268/A_2TB_Internal/RNA/RNA_lanes/Files_to_run_EBS/counts/Synon_counts/transpose.sh  | sed 's/0\\.5/3/g' > ", snp_dir, "merged_flt_c1-10.hmp321.RNAsetMAFover0under005_w_in_5kb_of_allgenes.KNNi.hmp.numeric_hetAs3.txt", sep="")

system(string1)
system(string3)
system(string5)
