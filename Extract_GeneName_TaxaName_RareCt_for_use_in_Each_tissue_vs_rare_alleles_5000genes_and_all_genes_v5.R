#Karl Kremling
#Feb 1 2017
#Extract number of rare upstream alleles per taxon per gene from a numericalized hmp file
#Run multiple instances on separate batches of genes (see below) of this script to speed it up, then combine batches at the end

library("miscTools")
library("readr")
library("reshape")
library("doParallel")
#install.packages("registerDoMc")
#library("registerDoMc")

#install.packages("doMC")
library("doMC")
registerDoMC(cores=48)

gff_file=read.delim("/local/workdir/kak268/rare_allel_calc/Zea_mays_gene_only.AGPv3.29.gff3", header=F, sep="\t", quote="") # add skip=6 if it's a full gff file

gff_file=read.delim("/home/kak268/work_data/maize_genome/Zea_mays_gene_only.AGPv3.29.gff3", header=F, sep="\t", quote="") # add skip=6 if it's a full gff file


colnames(gff_file)=c("Chr", "Annot_source", "Feature_type", "Start", "End", "blankcol1", "Strand", "blankcol2", "Feature_Name")
transcript_file=gff_file[gff_file$Feature_type=='gene',] # keep only the entries in the GFF file which are transcripts/genes
transcript_file$Feature_Name=sub(pattern= "ID=gene:", replacement='', x=transcript_file$Feature_Name)
transcript_file$Feature_Name=sub(pattern= ";.*", replacement='', x=transcript_file$Feature_Name)
transcript_file_only_chr1_to_chr10=transcript_file[transcript_file$Chr %in% c(1,2,3,4,5,6,7,8,9,10),] # keeps all transcript file info
transcript_file_only_chr1_to_chr10_gene_list=transcript_file[transcript_file$Chr %in% c(1,2,3,4,5,6,7,8,9,10),"Feature_Name"] # keeps only the genes on one of the ten chromosomes

#This was done in parallelat CBSU 
basedir="/local/workdir/kak268/rare_allel_calc/"
HMP_file=read.table(paste(basedir, "/merged_flt_c1-10.hmp321.RNAsetMAFover0under005_w_in_5kb_of_allgenes.KNNi.hmp.numeric_hetAs3.txt",sep=""), header=T)

rownames(HMP_file)=HMP_file[,1]
HMP_file=HMP_file[,-1]
HMP_file_w_pos=cbind(Chr=as.numeric(as.character(sub(pattern="S", replacement='', x=sub(pattern="_.*", replacement='', x=rownames(HMP_file))))), pos=as.numeric(as.character(sub(pattern=".*_", replacement='', x=rownames(HMP_file)))), as.data.frame(HMP_file))
#converts the list of SNPs to consistent format (otherwise they can be 1.0 or 1 etc and <NA>)
#SumNumChar=function(vec, val){sum(vec==val, na.rm=T)}

#divide HMP_file_tpose_w_pos into chr
for (i in 1:10){
  assign(paste("HMP_file_tpose_w_pos_Chr", i, sep=""), HMP_file_w_pos[HMP_file_w_pos$Chr==i,])  
}

#By breaking this into batches you can extract rare allele counts for multiple subsets of the genes at once in separate R sessions
batch=4
batch_points=c(0,4350,8700,13050,17400,21749,26099,30449,34799,39149)
registerDoMC(cores=12)

#count number of rare alleles upstream of each gene
parallel_5kb_rare_ct=NULL
for (g in c(0)){#, 1000, 2000, 3000, 4000)){
  prox_bound=g
  dist_bound=g+5000
  text_bound=paste(prox_bound/1000, "kb_to_", dist_bound/1000, "kb", sep="")
  cat("prox_bound is", prox_bound, "\n")
  cat("dist_bound is", dist_bound, "\n")
  cat("text_bound is", text_bound, "\n")
  parallel_5kb_rare_ct=NULL
  counter=0
  first_gene=batch_points[batch]+1
  sec_gene=batch_points[batch+1]
  
  for (i in unique(transcript_file_only_chr1_to_chr10$Feature_Name)[first_gene:sec_gene]){
    #i=unique(ex_bins_w_5kb_rare_ct$Feature_Name)[g]
    cat(".", sep="")
    if (counter%%50==0){
      cat(counter)
    }
    counter=counter+1
    #sub_melt=ex_bins_w_5kb_rare_ct[ex_bins_w_5kb_rare_ct$Feature_Name==i,]
    
    line=transcript_file_only_chr1_to_chr10[transcript_file_only_chr1_to_chr10$Feature_Name==i,]
    SNPs=NULL 
    #undo#GERP_subset=NULL
    #GERP_subset=min_2_GERPs[min_2_GERPs$Chr==line$Chr,]
    HMP_file_tpose_w_pos_single_Chr=get(paste("HMP_file_tpose_w_pos_Chr", line$Chr, sep=""))
    #undo#min_2_GERPs_single_Chr=get(paste("min_2_GERPs_Chr", line$Chr, sep=""))
    if (line$Strand=="-"){
      fivePend=line$End
      SNPs=HMP_file_tpose_w_pos_single_Chr[HMP_file_tpose_w_pos_single_Chr$pos<=fivePend+dist_bound & HMP_file_tpose_w_pos_single_Chr$pos>fivePend+prox_bound,] # getting the SNPs 5kb upstream
      #undo#GERP_subset=min_2_GERPs_single_Chr[min_2_GERPs_single_Chr$pos<=fivePend+5000 & min_2_GERPs_single_Chr$pos>fivePend,]
    } else if (line$Strand=="+"){
      fivePend=line$Start
      SNPs=HMP_file_tpose_w_pos_single_Chr[HMP_file_tpose_w_pos_single_Chr$pos>=fivePend-dist_bound & HMP_file_tpose_w_pos_single_Chr$pos<fivePend-prox_bound,]
      #undo#GERP_subset=min_2_GERPs_single_Chr[min_2_GERPs_single_Chr$pos>=fivePend-5000 & min_2_GERPs_single_Chr$pos<fivePend,]
    }
    temp_parallel_5kb_rare_ct=data.frame(Taxon=NA, Feature_Name=NA, Upstr_5kb_rare_homo_ct=NA)#, Upstr_5kb_rare_het_ct=NA)
    
    temp_parallel_5kb_rare_ct=foreach(k=3:(length(colnames(SNPs))), .combine=rbind) %dopar% { 
      #for(j in sub_melt$Taxon){
      j=colnames(SNPs)[k] # only if using foreach in this loop
      #single_taxon_SNPs=as.numeric(as.character(SNPs[,j]))
      #cat(j)
      #orig#single_taxon_SNPs=as.data.frame(cbind(Chr=as.numeric(as.character(SNPs$Chr)), pos=as.numeric(as.character(SNPs$pos)), Geno=as.numeric(as.character(SNPs[,j]))))
      single_taxon_SNPs=as.data.frame(cbind(pos=SNPs$pos, Geno=SNPs[,j]))
      #ex_bins_w_5kb_rare_ct[ex_bins_w_5kb_rare_ct$Feature_Name==i & ex_bins_w_5kb_rare_ct$Taxon==j,"Upstr_5kb_rare_homo_ct"]=sum(single_taxon_SNPs$Geno==0, na.rm=T)
      #ex_bins_w_5kb_rare_ct[ex_bins_w_5kb_rare_ct$Feature_Name==i & ex_bins_w_5kb_rare_ct$Taxon==j,"Upstr_5kb_rare_het_ct"]=sum(single_taxon_SNPs$Geno==3, na.rm=T)
      #Record number of rare SNPs in GERP>2 site
      #undo#single_taxon_SNPs_w_GERP=merge(single_taxon_SNPs, GERP_subset, by="pos")
      #undo#ex_bins_w_5kb_rare_ct[ex_bins_w_5kb_rare_ct$Feature_Name==i & ex_bins_w_5kb_rare_ct$Taxon==j,c("Upstr_5kb_rare_homo_ct", "Upstr_5kb_rare_het_ct", "Upstr_5kb_rare_homo_ct_GERP_over_2", "Upstr_5kb_rare_het_ct_GERP_over_2")]=c(sum(single_taxon_SNPs$Geno==0, na.rm=T), sum(single_taxon_SNPs$Geno==3, na.rm=T), sum(single_taxon_SNPs_w_GERP$Geno==0, na.rm=T), sum(single_taxon_SNPs_w_GERP$Geno==3, na.rm=T))
      data.frame(i, j, sum(single_taxon_SNPs$Geno==0, na.rm=T))#, sum(single_taxon_SNPs$Geno==3, na.rm=T))
      
      #ex_bins_w_5kb_rare_ct[ex_bins_w_5kb_rare_ct$Feature_Name==i & ex_bins_w_5kb_rare_ct$Taxon==j,"Upstr_5kb_rare_homo_ct_GERP_over_2"]=sum(single_taxon_SNPs_w_GERP$Geno==0, na.rm=T)
      #ex_bins_w_5kb_rare_ct[ex_bins_w_5kb_rare_ct$Feature_Name==i & ex_bins_w_5kb_rare_ct$Taxon==j,"Upstr_5kb_rare_het_ct_GERP_over_2"]=sum(single_taxon_SNPs_w_GERP$Geno==3, na.rm=T)
    }
    if (is.null(parallel_5kb_rare_ct)){
      parallel_5kb_rare_ct=temp_parallel_5kb_rare_ct
    } else {
      parallel_5kb_rare_ct=rbind(parallel_5kb_rare_ct,temp_parallel_5kb_rare_ct)
    }
  }
  date()
  colnames(parallel_5kb_rare_ct)=c("Feature_Name","Taxon", "Upstr_5kb_rare_homo_ct")#, "Upstr_5kb_rare_het_ct")
  
  #ex_bins_w_5kb_rare_ct=merge(ex_bins_w_5kb_rare_ct, parallel_5kb_rare_ct, by=c("Feature_Name","Taxon"))
  parallel_5kb_rare_ct$Taxon=sub(pattern="X282set", replacement="282set", parallel_5kb_rare_ct$Taxon)
  parallel_5kb_rare_ct$Taxon=gsub(pattern="\\.", replacement="-", parallel_5kb_rare_ct$Taxon)
  cat("head of ex_bins_w_5kb_rare_ct is :")
  head(parallel_5kb_rare_ct)
  write.table(parallel_5kb_rare_ct, file=paste(basedir,"5kb_rare_ct_all_genes_w_hmp321_taxa_",batch,".txt", sep=""), quote=F, row.names = F)
}


#join separate rare allele file count subsets
full_file=NULL
basedir="/local/workdir/kak268/rare_allel_calc/"
for (batch in 1:8){
  cat("batch: ", batch, "\n")
  current_file=read.table(file=paste(basedir,"5kb_rare_ct_all_genes_w_hmp321_taxa_",batch,".txt", sep=""), header = T)
  cat("current_file dim is ", dim(current_file),"\n")
  
  #if(exists("full_file")){
    full_file=rbind(full_file, current_file)
  #}else{
  #  cat(" in else \n")
  #  full_file=current_file
  #}
}
write.table(full_file, file="5kb_rare_ct_all_genes_w_hmp321_taxa_all_parts.txt", quote=F, row.names = F)

#convert rare allele counts to a 1 if there is >=1 rare allele present, but 0 if not
#this is the binary "rare allele present upstream" vs "no rare allele upstream" test done as part of the response to reviewers
full_file_pres_abs=full_file

full_file_pres_abs[full_file_pres_abs$Upstr_5kb_rare_homo_ct>=1,3]=1

write.table(full_file_pres_abs, file=paste(basedir,"5kb_rare_pres_abs_all_genes_w_hmp321_taxa_all_parts.txt", sep=""), quote=F, row.names = F)



