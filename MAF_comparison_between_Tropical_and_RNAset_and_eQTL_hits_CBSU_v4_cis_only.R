#!/bin/Rscript
#CAll list this: Rscript --verbose MAF_comparison_between_Tropical_and_RNAset_and_eQTL_hits_CBSU.R GShoot
#or make bash script as follows and run with perl_fork_univ.pl
#!/usr/bin/env bash   
#Rscript MAF_comparison_between_Tropical_and_RNAset_and_eQTL_hits_CBSU_v4_cis_only.R L3Tip
#Rscript MAF_comparison_between_Tropical_and_RNAset_and_eQTL_hits_CBSU_v4_cis_only.R L3Base
#Rscript MAF_comparison_between_Tropical_and_RNAset_and_eQTL_hits_CBSU_v4_cis_only.R GShoot
#Rscript MAF_comparison_between_Tropical_and_RNAset_and_eQTL_hits_CBSU_v4_cis_only.R GRoot
#Rscript MAF_comparison_between_Tropical_and_RNAset_and_eQTL_hits_CBSU_v4_cis_only.R LMAD
#Rscript MAF_comparison_between_Tropical_and_RNAset_and_eQTL_hits_CBSU_v4_cis_only.R LMAN
#Rscript MAF_comparison_between_Tropical_and_RNAset_and_eQTL_hits_CBSU_v4_cis_only.R Kern

#perl /programs/bin/perlscripts/perl_fork_univ.pl masterv4.txt 7
#mkdir /SSD/kak268
#mkdir /SSD/kak268/genotypes/
#mkdir /SSD/kak268/genotypes/MAF_summaries/
#mkdir /SSD/kak268/genotypes/plots/
# mkdir /SSD/kak268/genotypes/top_eQTL_results/
cbind.fill <- function(...){
  nm <- list(...) 
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

args <- commandArgs(TRUE)
#CorrAns = args[1]
i = as.character(args[1])


basedir="/SSD/kak268/genotypes/"
eqtldir="/SSD/kak268/genotypes/top_eQTL_results/"
plotdir=paste("/SSD/kak268/genotypes/plots_v5/", i,"/", sep = "")
dir.create(plotdir)
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
#CorrAns = args[1]
#i="GRoot"
#for (i in c("GRoot", "Kern", "L3Base", "L3Tip", "GShoot", "LMAD", "LMAN")){
library("ggplot2")
library("plyr")
library("readr")
library("RColorBrewer")
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
#multiplot http://rstudio-pubs-static.s3.amazonaws.com/2852_379274d7c5734f979e106dcf019ec46c.html
#multiplot http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/

#Trop_MAFs=read.table("/workdir/kak268/genotypes/allchr_summary3Tropical_merged_flt_c1-10.hmp321.KNNi.hmp.txt", header=T)
Trop_MAFs=read_delim(paste(basedir,"MAF_summaries/allchr_summary3Tropical_merged_flt_c1-10.hmp321.KNNi.hmp.txt", sep=""), col_names=T, delim="\t")
#5 SNPS are duplicated because they exist in multiple staes in the hmp geno summary file, thus I am excluding them: 2 S8_173883313 2 S6_165353519 2 S3_226897037 2 S2_233507851 2 S10_147291022 
Trop_MAFs=Trop_MAFs[!duplicated(Trop_MAFs$SNPID),]
#for ( i in c("L3Tip", "Kern", "GRoot", "GShoot", "LMAD", "LMAN", "L3Base")){
#i="LMAN"
cat(i)
cat(date(), "\n")
#Read in the MAFs
#TopTrop_MAFs=Trop_MAFs[1:1000000,]
###Read in the RNAset material MAFs
#RNASet_MAFs=read.table(paste("/workdir/kak268/genotypes/MAF_summaries/allchr_site_summary_",i,"_mdp_genotype.txt", sep = ""), sep = "\t", header=T)
RNASet_MAFs=read_delim(paste(basedir, "MAF_summaries/allchr_site_summary_",i,"_mdp_genotype.txt", sep = ""), delim = "\t", col_names = T)
cat("dim of RNASet_MAFs is ", dim(RNASet_MAFs), " ")
#ensure you don't plot any dots with RNAset maf under 0.05
RNASet_MAFs=RNASet_MAFs[which(RNASet_MAFs$MAF>=0.05),]
#}#end of dim loop if only trying to count number of SNPs in RNAsets and 

#TopRNASet_MAFs=RNASet_MAFs[1:1000000,]
#merged_MAFs=merge(Trop_MAFs, RNASet_MAFs, by="SNPID")
merged_MAFs=join(Trop_MAFs, RNASet_MAFs, by="SNPID", type="inner")
#merged_MAFs_no_first=join(Trop_MAFs, RNASet_MAFs, by="SNPID", type="inner")
#merged_MAFs_no_first_or_inner=join(Trop_MAFs, RNASet_MAFs, by="SNPID")

colnames(merged_MAFs)=c("SNPID", "Tropical_MAF", "RNASet_MAF")
cat("dim of merged_MAFs is ", dim(merged_MAFs), " ")
cat("pct over or equal 0.2 RNAset, and tropical under 0.05 ", dim(merged_MAFs[which(merged_MAFs$Tropical_MAF<0.05 & merged_MAFs$RNASet_MAF>=0.2),]))


#version of eqtl hit file before avg'ing dups
#eQTL_results=read.table(paste(eqtldir,"top_hit_by_SNP_", i, "_merged_flt_allchr_maxP000001.hmp321.onlyRNAset_MAFover005.KNNi.w.df_STAR_HTSeq_counts_B73_match_based_on_genet_dist_DESeq2_normed_rounded_origNames_and_Altnames_sm_rand_and_box_coxed_w_MDSPCs_and_25PEERs.txt", sep = ""), header = T)
#eqtl hits with dups avged
#eQTL_results=read.table(paste(eqtldir,"top_hit_by_SNP_", i, "_merged_flt_allchr_maxP000001.hmp321.onlyRNAset_MAFover005.KNNi.w.df_STAR_HTSeq_counts_B73_match_based_on_genet_dist_DESeq2_normed_rounded_origNames_and_Altnames_sm_rand_and_box_coxed_w_MDSPCs_and_25PEERs_avgDups.txt", sep = ""), header = T)
eQTL_results=read_delim(paste(eqtldir,"top_hit_by_SNP_", i, "_merged_flt_allchr_maxP000001.hmp321.onlyRNAset_MAFover005.KNNi.w.df_STAR_HTSeq_counts_B73_match_based_on_genet_dist_DESeq2_normed_rounded_origNames_and_Altnames_sm_rand_and_box_coxed_w_MDSPCs_and_25PEERs_avgDups.txt", sep = ""), delim = "\t", col_names = T)
colnames(eQTL_results)[2]="SNPID"
## cant do this due to abnormal characters eQTL_results=read_delim(paste("/workdir/kak268/eQTL_results/top_eQTL_results/top_hit_by_SNP_", i, "_merged_flt_allchr_maxP000001.hmp321.onlyRNASet_MAFover005.KNNi.w.df_STAR_HTSeq_counts_B73_match_based_on_genet_dist_DESeq2_normed_rounded_origNames_and_Altnames_sm_rand_and_box_coxed_w_MDSPCs_and_25PEERs.txt", sep = ""), delim = "\t", col_names = T)
cat("dim of eQTL_results is ", dim(eQTL_results), " ")

merged_MAFs_w_eQTL=join(merged_MAFs, eQTL_results, by = "SNPID", type="inner")
#merged_MAFs_w_eQTL_no_first=join(merged_MAFs, eQTL_results, by = "SNPID", type="inner")
#merged_MAFs_w_eQTL_no_first_no_inner=join(merged_MAFs, eQTL_results, by = "SNPID")

merged_MAFs_w_eQTL$r2=as.numeric(as.character(merged_MAFs_w_eQTL$r2))
cat("dim of merged_MAFs_w_eQTL is ", dim(merged_MAFs_w_eQTL), " ")
merged_MAFs_w_eQTL_sampled=merged_MAFs_w_eQTL[sample(nrow(merged_MAFs_w_eQTL),size=3000000,replace=FALSE),]


merged_MAFs_w_eQTL$Trop_MAF_bins= round(merged_MAFs_w_eQTL$Tropical_MAF,digits = 2)
bin_Count=as.data.frame(unclass(rle(x = sort(merged_MAFs_w_eQTL$Trop_MAF_bins)))) # calcualtes the number of observations in each bin
merged_MAFs_w_eQTL_count = merge(merged_MAFs_w_eQTL, bin_Count, by.x="Trop_MAF_bins", by.y="values")
#merged_MAFs_w_eQTL_count$lengths=1 # include this option if you don't want to color by number of SNPs per bin

# #Start CHECK
# 
#1D Trop MAF plot with R2 from RNASet MAF>0.2 colored by bin count
tiff(paste(plotdir,i,"_v_R2_in_RNASetMAF_v_", i, "_top_eQTL_per_SNP_allchr.tiff", sep=""), height=800, width=2200, res=800)
MAF_over_0.2_merged_MAFs_w_eQTL_count=subset(merged_MAFs_w_eQTL_count, RNASet_MAF>=0.2)
#write.table(x=MAF_over_0.2_merged_MAFs_w_eQTL_count, file = paste(plotdir,i,"_v_R2_in_RNASetMAF_v_", i, "_top_eQTL_per_SNP_allchr.txt", sep=""), row.names = F)
plotbar1=ggplot(data=MAF_over_0.2_merged_MAFs_w_eQTL_count, aes(x = Trop_MAF_bins, y = r2, fill=lengths), stat="identity") +
  stat_summary(fun.y = "mean", geom = "bar") +
  scale_fill_gradientn(colours=colorspace::diverge_hcl(n = 8, l=c(30,90), c=c(80,30))) +
  #scale_fill_gradientn(colours=colorspace::heat_hcl(n = 8, l=c(30,90))) +
  #scale_fill_gradientn(colours=colorspace::rainbow_hcl(n = 8, l=c(30,90))) +
  theme(axis.text=element_text(size=2), axis.title=element_text(size=4), legend.title=element_text(size=4), legend.text=element_text(size=2), legend.key.size=unit(0.15,"cm"),
        panel.grid.major = element_blank(), axis.title.y=element_text(margin=margin(0,5,0,0)), axis.title.x=element_text(margin=margin(5,0,0,0)), panel.grid.minor = element_blank(), panel.border = element_blank(), axis.ticks.length = unit(0.1, "cm"), axis.ticks=element_line(color="black"), panel.background = element_blank(), axis.line = element_line(colour = "black", size = 0.3))  +
  labs(x = "Tropical MAF bins", y = expression(RNA~ set~ eQTL~ mean~ R^2), fill="SNPs per tropical \nMAF bin")
print(plotbar1)
dev.off()
# 
# 
# #1D Trop MAF plot with R2 from RNASet MAF>0.2 NOT colored by bin count
# tiff(paste(plotdir, "OneD_one_color_TropMAF_v_",i,"v_R2_in_RNASetMAF_v_", i, "_top_eQTL_per_SNP_allchr.tiff", sep=""), height=800, width=2200, res=800)
# plotbar2=ggplot(data=subset(merged_MAFs_w_eQTL_count, RNASet_MAF>=0.2), aes(x = Trop_MAF_bins, y = r2, fill=lengths), stat="identity", color="dodgerblue") +
#   stat_summary(fun.y = "mean", geom = "bar") +
#   #scale_fill_gradientn(colours=colorspace::diverge_hcl(n = 8, l=c(30,90), c=c(80,30))) +
#   #scale_fill_gradientn(colours=colorspace::heat_hcl(n = 8, l=c(30,90))) +
#   #scale_fill_gradientn(colours=colorspace::rainbow_hcl(n = 8, l=c(30,90))) +
#   theme(axis.text=element_text(size=2), axis.title=element_text(size=4), legend.title=element_text(size=4), legend.text=element_text(size=2), legend.key.size=unit(0.15,"cm"),
#         panel.grid.major = element_blank(), axis.title.y=element_text(margin=margin(0,5,0,0)), axis.title.x=element_text(margin=margin(5,0,0,0)), panel.grid.minor = element_blank(), panel.border = element_blank(), axis.ticks.length = unit(0.1, "cm"), axis.ticks=element_line(color="black"), panel.background = element_blank(), axis.line = element_line(colour = "black", size = 1))  +
#   labs(x = "Tropical MAF bins", y = expression(RNA~ set~ eQTL~ mean~ R^2), fill="SNPs per tropical \nMAF bin")
# print(plotbar2)
# dev.off()
# #Stop CHECK
# 
# # Start Check
#color NOT scaled by R2 2D MAF plot
# tiff(paste(plotdir,"TwoD_one_color_TropMAF_v_",i,"RNASetMAF_v_", i, "_top_eQTL_per_SNP_allchr.tiff", sep=""), height=1688, width=2700, res=500)
# #plot=qplot(aes_string(x="Tropical_MAF", y=paste(i,"_RNASet_MAFs", sep="")), data = merged_MAFs_w_eQTL, colour=r2, alpha=I(0.015), size=I(0.5)) +
# #plot=qplot(x=Tropical_MAF, y=eval(parse(text=paste(i,"_RNASet_MAFs", sep=""))), data = merged_MAFs_w_eQTL, colour=r2, alpha=I(0.025), size=I(1)) +
# #plot=qplot(x=Tropical_MAF, y=paste(i,"_RNASet_MAFs", sep=""), data = merged_MAFs_w_eQTL, colour=r2, alpha=I(0.025), size=I(1)) +
# plot1=qplot(Tropical_MAF, RNASet_MAF, data = merged_MAFs_w_eQTL_sampled, colour=r2, alpha=I(0.02), size=I(0.01)) +
#   scale_color_continuous(low = "dodgerblue", high="dodgerblue", name=bquote(R^{2})) +
#   #coord_fixed() +
#   coord_fixed(xlim = c(0, 0.5), ylim = c(0, 0.5), ratio=1) +
#   theme(axis.text=element_text(size=7), axis.title=element_text(size=9), legend.title=element_text(size=9), legend.text=element_text(size=9), legend.key.size=unit(0.4,"cm"), panel.grid.major = element_blank(), axis.title.y=element_text(margin=margin(0,5,0,0)), axis.title.x=element_text(margin=margin(5,0,0,0)), panel.grid.minor = element_blank(), panel.border = element_blank(), axis.ticks.length = unit(0.1, "cm"), axis.ticks=element_line(color="black"), panel.background = element_blank(), axis.line = element_line(colour = "black", size = 1)) 
# print(plot1)
# dev.off()


# #START Check
##color scaled by R2 2D MAF plots
tiff(paste(plotdir,"TwoD_TropMAF_v_",i,"RNASetMAF_v_", i, "_top_eQTL_per_SNP_allchr.tiff", sep=""), height=1688, width=2700, res=500)
#plot=qplot(aes_string(x="Tropical_MAF", y=paste(i,"_RNASet_MAFs", sep="")), data = merged_MAFs_w_eQTL, colour=r2, alpha=I(0.015), size=I(0.5)) +
#plot=qplot(x=Tropical_MAF, y=eval(parse(text=paste(i,"_RNASet_MAFs", sep=""))), data = merged_MAFs_w_eQTL, colour=r2, alpha=I(0.025), size=I(1)) +
#plot=qplot(x=Tropical_MAF, y=paste(i,"_RNASet_MAFs", sep=""), data = merged_MAFs_w_eQTL, colour=r2, alpha=I(0.025), size=I(1)) +
plot2=qplot(Tropical_MAF, RNASet_MAF, data = merged_MAFs_w_eQTL_sampled, colour=r2, alpha=I(0.02), size=I(0.01)) +
  scale_color_continuous(low = "dodgerblue", high="violetred3", name=bquote(R^{2})) + # violetred3
  coord_fixed(xlim = c(0, 0.5), ylim = c(0, 0.5), ratio=1) +
  #coord_fixed() +
  theme(axis.text=element_text(size=7), axis.title=element_text(size=9), legend.title=element_text( size=9), legend.text=element_text(size=9), legend.key.size=unit(0.4,"cm"),
        panel.grid.major = element_blank(), axis.title.y=element_text(margin=margin(0,5,0,0)), axis.title.x=element_text(margin=margin(5,0,0,0)), panel.grid.minor = element_blank(), panel.border = element_blank(), axis.ticks.length = unit(0.1, "cm"), axis.ticks=element_line(color="black"), panel.background = element_blank(), axis.line = element_line(colour = "black", size = 0.3))
print(plot2)
dev.off()
# Stop Check

#START CHECK
#print side by side histogram distributions
for (upper_bin_min in as.numeric(c(0, 0.1, 0.2, 0.3, 0.4))){
  cat(upper_bin_min, " ", is.numeric(upper_bin_min))
  #on left edge
  data1=subset(merged_MAFs_w_eQTL, RNASet_MAF>=upper_bin_min & RNASet_MAF<upper_bin_min+0.1 & Tropical_MAF>=0 & Tropical_MAF<0.1)
  #on diag
  data2=subset(merged_MAFs_w_eQTL, RNASet_MAF>=upper_bin_min & RNASet_MAF<upper_bin_min+0.1 & Tropical_MAF>=upper_bin_min & Tropical_MAF<upper_bin_min+0.1)
  data3=cbind.fill(data1, data2)
  #write.table(data3, file = paste(plotdir,"Hist_0_and_", upper_bin_min, "_bins_TropMAF_v_RNASetMAF_v_", i, "_top_eQTL_per_SNP.txt", sep=""), row.names = F)
  tiff(compression = "lzw", filename = paste(plotdir,"Hist_0_and_", upper_bin_min, "_bins_TropMAF_v_RNASetMAF_v_", i, "_top_eQTL_per_SNP.tiff", sep=""), height=1688, width=2700, res=800)
  histgg1=ggplot(merged_MAFs_w_eQTL, aes(r2)) +
    geom_histogram(aes(y=..density..),data1, fill="darkred", alpha=0.5) +
    geom_histogram(aes(y=..density..),data2, fill="midnightblue", alpha=0.5) +
    theme(axis.text=element_text(size=7), axis.title=element_text(size=9), legend.title=element_text(size=9), legend.text=element_text(size=9), legend.key.size=unit(0.4,"cm"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black", size = 0.3)) +
    xlab(bquote(R^{2}))+
    ggtitle(paste("Wilc. s. rank P = ", signif(wilcox.test(data1$r2, data2$r2)$p.value, 6), "\n KS P =", signif(ks.test(data1$r2, data2$r2)$p.value, 6), sep=""))
  print(histgg1)
  dev.off()
}
# # Stop CHECK

plots=list()
for (chr in c(1,2,3,4,5,6,7,8,9,10)){
  cat("chr is ", chr)
  
  # #Start Check
  #1D Trop MAF plot with R2 from RNAset MAF>0.2 Colored by counts per bin
  one_chr_merged_MAFs_w_eQTL=subset(merged_MAFs_w_eQTL, Chr==chr)
  tiff(paste(plotdir,"OneD_chr_",chr,"_TropMAF_v_",i,"v_R2_in_RNASetMAF_v_", i, "_top_eQTL_per_SNP_allchr.tiff", sep=""), height=800, width=2200, res=800)
  one_chr_merged_MAFs_w_eQTL$Trop_MAF_bins= round(one_chr_merged_MAFs_w_eQTL$Tropical_MAF,digits = 2)
  bin_Count=as.data.frame(unclass(rle(x = sort(one_chr_merged_MAFs_w_eQTL$Trop_MAF_bins)))) # calcualtes the number of observations in each bin
  one_chr_merged_MAFs_w_eQTL_count = merge(one_chr_merged_MAFs_w_eQTL, bin_Count, by.x="Trop_MAF_bins", by.y="values")
  MAF_over_0.2_one_chr_merged_MAFs_w_eQTL_count=subset(one_chr_merged_MAFs_w_eQTL_count, RNASet_MAF>=0.2)
  #write.table(x=MAF_over_0.2_one_chr_merged_MAFs_w_eQTL_count, file = paste(plotdir,"OneD_chr_",chr,"_TropMAF_v_",i,"v_R2_in_RNASetMAF_v_", i, "_top_eQTL_per_SNP_allchr.tiff", sep=""), row.names = F, quote = F)
  plotbar3=ggplot(data=MAF_over_0.2_one_chr_merged_MAFs_w_eQTL_count, aes(x = Trop_MAF_bins, y = r2, fill=lengths), stat="identity") +
    stat_summary(fun.y = "mean", geom = "bar") +
    scale_fill_gradientn(colours=colorspace::diverge_hcl(n = 8, l=c(30,90), c=c(80,30))) +
    #scale_fill_gradientn(colours=colorspace::heat_hcl(n = 8, l=c(30,90))) +
    #scale_fill_gradientn(colours=colorspace::rainbow_hcl(n = 8, l=c(30,90))) +
    theme(axis.text=element_text(size=2), axis.title=element_text(size=4), legend.title=element_text(size=4), legend.text=element_text(size=2), legend.key.size=unit(0.15,"cm"),
          panel.grid.major = element_blank(), axis.title.y=element_text(margin=margin(0,5,0,0)), axis.title.x=element_text(margin=margin(5,0,0,0)), panel.grid.minor = element_blank(), panel.border = element_blank(), axis.ticks.length = unit(0.1, "cm"), axis.ticks=element_line(color="black"), panel.background = element_blank(), axis.line = element_line(colour = "black", size = 0.3))  +
    labs(x = "Tropical MAF bins", y = expression(RNA~ set~ eQTL~ mean~ R^2), fill="SNPs per tropical \nMAF bin")
  print(plotbar3)
  dev.off()
  # #1D Trop MAF plot with R2 from RNAset MAF>0.2 NOT Colored by counts per bin
  # one_chr_merged_MAFs_w_eQTL=subset(merged_MAFs_w_eQTL, Chr==chr)
  # tiff(paste(plotdir,"OneD_one_color_chr_",chr,"_TropMAF_v_",i,"v_R2_in_RNASetMAF_v_", i, "_top_eQTL_per_SNP_allchr.tiff", sep=""), height=800, width=2200, res=800)
  # one_chr_merged_MAFs_w_eQTL$Trop_MAF_bins= round(one_chr_merged_MAFs_w_eQTL$Tropical_MAF,digits = 2)
  # bin_Count=as.data.frame(unclass(rle(x = sort(one_chr_merged_MAFs_w_eQTL$Trop_MAF_bins)))) # calcualtes the number of observations in each bin
  # one_chr_merged_MAFs_w_eQTL_count = merge(one_chr_merged_MAFs_w_eQTL, bin_Count, by.x="Trop_MAF_bins", by.y="values")
  # one_chr_merged_MAFs_w_eQTL_count$lengths=1
  # plotbar4=ggplot(data=subset(one_chr_merged_MAFs_w_eQTL_count, RNASet_MAF>=0.2), aes(x = Trop_MAF_bins, y = r2, fill=lengths), stat="identity", color="dodgerblue") +
  #   stat_summary(fun.y = "mean", geom = "bar") +
  #   #scale_fill_gradientn(colours=colorspace::diverge_hcl(n = 8, l=c(30,90), c=c(80,30))) +
  #   #scale_fill_gradientn(colours=colorspace::heat_hcl(n = 8, l=c(30,90))) +
  #   #scale_fill_gradientn(colours=colorspace::rainbow_hcl(n = 8, l=c(30,90))) +
  #   theme(axis.text=element_text(size=2), axis.title=element_text(size=4), legend.title=element_text(size=4), legend.text=element_text(size=2), legend.key.size=unit(0.15,"cm"),
  #         panel.grid.major = element_blank(), axis.title.y=element_text(margin=margin(0,5,0,0)), axis.title.x=element_text(margin=margin(5,0,0,0)), panel.grid.minor = element_blank(), panel.border = element_blank(), axis.ticks.length = unit(0.1, "cm"), axis.ticks=element_line(color="black"), panel.background = element_blank(), axis.line = element_line(colour = "black", size = 1))  +
  #   labs(x = "Tropical MAF bins", y = expression(RNA~ set~ eQTL~ mean~ R^2), fill="SNPs per tropical \nMAF bin")
  # print(plotbar4)
  # dev.off()
  # # Stop CHeck
  
  # # Start Check
  #color not scaled by R2 2D MAF plot
  # tiff(paste(plotdir,"TwoD_chr_",chr,"_one_color_TropMAF_v_",i,"RNASetMAF_v_", i, "_top_eQTL_per_SNP_allchr.tiff", sep=""), height=1688, width=2700, res=500)
  # #plot=qplot(aes_string(x="Tropical_MAF", y=paste(i,"_RNASet_MAFs", sep="")), data = merged_MAFs_w_eQTL, colour=r2, alpha=I(0.015), size=I(0.5)) +
  # #plot=qplot(x=Tropical_MAF, y=eval(parse(text=paste(i,"_RNASet_MAFs", sep=""))), data = merged_MAFs_w_eQTL, colour=r2, alpha=I(0.025), size=I(1)) +
  # #plot=qplot(x=Tropical_MAF, y=paste(i,"_RNASet_MAFs", sep=""), data = merged_MAFs_w_eQTL, colour=r2, alpha=I(0.025), size=I(1)) +
  # plot3=qplot(Tropical_MAF, RNASet_MAF, data = subset(merged_MAFs_w_eQTL, Chr==chr), colour=r2, alpha=I(0.02), size=I(0.01)) +
  #   scale_color_continuous(low = "dodgerblue", high="dodgerblue", name=bquote(R^{2})) +
  #   #coord_fixed() +
  #   coord_fixed(xlim = c(0, 0.5), ylim = c(0, 0.5), ratio=1) +
  #   theme(axis.text=element_text(size=7), axis.title=element_text(size=9), legend.title=element_text( size=9), legend.text=element_text(size=9), legend.key.size=unit(0.4,"cm"),
  #         panel.grid.major = element_blank(), axis.title.y=element_text(margin=margin(0,5,0,0)), axis.title.x=element_text(margin=margin(5,0,0,0)), panel.grid.minor = element_blank(), panel.border = element_blank(), axis.ticks.length = unit(0.1, "cm"), axis.ticks=element_line(color="black"), panel.background = element_blank(), axis.line = element_line(colour = "black", size = 1)) 
  # print(plot3)
  # dev.off()
  # # Stop Check
  
  
  #color scaled by R2  2D MAF plot
  ##CHECK tiff(paste(plotdir,"TwoD_chr_",chr,"_TropMAF_v_",i,"RNASetMAF_v_", i, "_top_eQTL_per_SNP_allchr.tiff", sep=""), height=1688, width=2700, res=500)
  #plot=qplot(aes_string(x="Tropical_MAF", y=paste(i,"_RNASet_MAFs", sep="")), data = merged_MAFs_w_eQTL, colour=r2, alpha=I(0.015), size=I(0.5)) +
  #plot=qplot(x=Tropical_MAF, y=eval(parse(text=paste(i,"_RNASet_MAFs", sep=""))), data = merged_MAFs_w_eQTL, colour=r2, alpha=I(0.025), size=I(1)) +
  #plot=qplot(x=Tropical_MAF, y=paste(i,"_RNASet_MAFs", sep=""), data = merged_MAFs_w_eQTL, colour=r2, alpha=I(0.025), size=I(1)) +
  plot4=qplot(Tropical_MAF, RNASet_MAF, data = subset(merged_MAFs_w_eQTL, Chr==chr), colour=r2, alpha=I(0.02), size=I(0.01)) +
    scale_color_continuous(low = "dodgerblue", high="violetred3", name=bquote(R^{2})) + # violetred3
    coord_fixed(xlim = c(0, 0.5), ylim = c(0, 0.5), ratio=1) +
    #coord_fixed() +
    ggtitle(paste("Chr. ", chr," ", i ))+
    theme(axis.text=element_text(size=7), axis.title=element_text(size=9), legend.title=element_text( size=9), legend.text=element_text(size=9), legend.key.size=unit(0.4,"cm"),
          panel.grid.major = element_blank(), axis.title.y=element_text(margin=margin(0,5,0,0)), axis.title.x=element_text(margin=margin(5,0,0,0)), panel.grid.minor = element_blank(), panel.border = element_blank(), axis.ticks.length = unit(0.1, "cm"), axis.ticks=element_line(color="black"), panel.background = element_blank(), axis.line = element_line(colour = "black", size = 0.3)) 
  plots[[chr]] = plot4
  ##CHECK print(plot4)
  ##CHECK dev.off()
  
  
  #Start Check
  #Make chromosome based histograms of R2
  for (upper_bin_min in as.numeric(c(0, 0.1, 0.2, 0.3, 0.4))){
    cat(upper_bin_min, " ", is.numeric(upper_bin_min))
    #on left edge
    data1=subset(merged_MAFs_w_eQTL, Chr==chr & RNASet_MAF>=upper_bin_min & RNASet_MAF<upper_bin_min+0.1 & Tropical_MAF>=0 & Tropical_MAF<0.1)
    #on diag
    data2=subset(merged_MAFs_w_eQTL, Chr==chr & RNASet_MAF>=upper_bin_min & RNASet_MAF<upper_bin_min+0.1 & Tropical_MAF>=upper_bin_min & Tropical_MAF<upper_bin_min+0.1)
    data3=cbind.fill(data1, data2)
    #write.table(data3, file = paste(plotdir,"Hist_chr_",chr,"_0_and_", upper_bin_min, "_bins_TropMAF_v_RNASetMAF_v_", i, "_top_eQTL_per_SNP.txt", sep=""), quote=F, row.names = F)
    tiff(compression = "lzw", filename = paste(plotdir,"Hist_chr_",chr,"_0_and_", upper_bin_min, "_bins_TropMAF_v_RNASetMAF_v_", i, "_top_eQTL_per_SNP.tiff", sep=""), height=1688, width=2700, res=800)
    histgg5=ggplot(merged_MAFs_w_eQTL, aes(r2)) +
      geom_histogram(aes(y=..density..),data1, fill="darkred", alpha=0.5) +
      geom_histogram(aes(y=..density..),data2, fill="midnightblue", alpha=0.5) +
      theme(axis.text=element_text(size=7), axis.title=element_text(size=9), legend.title=element_text(size=9), legend.text=element_text(size=9), legend.key.size=unit(0.4,"cm"),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black", size = 0.3)) +
      xlab(bquote(R^{2}))+
      ggtitle(paste("Wilc. s. rank P = ", signif(wilcox.test(data1$r2, data2$r2)$p.value, 6), "\n KS P =", signif(ks.test(data1$r2, data2$r2)$p.value, 6), sep=""))
    print(histgg5)
    dev.off()
  }
  #Stop CHECK
}

#tiff(paste(plotdir,"TwoD_all_sep_chr_TropMAF_v_",i,"RNASetMAF_v_", i, "_top_eQTL_per_SNP_allchr.tiff", sep=""), height=8440, width=5400, res=500)
tiff(paste(plotdir,"TwoD_all_sep_chr_TropMAF_v_",i,"RNASetMAF_v_", i, "_top_eQTL_per_SNP_allchr.tiff", sep=""), height=8440, width=5400, res=400)
multiplot(plotlist=plots, cols=2)
dev.off()