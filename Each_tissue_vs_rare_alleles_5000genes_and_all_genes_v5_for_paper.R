#Karl Kremling
#Feb 2017
# Correlate expression rank with rare allele abundance for different subsets of expressed genes (ie top 5k next 5k, or 0-1k, 1-2k, 2-3k, etc) 


##
##
##Plot smile plots for all tissues with different fractions of top 10k subset, also has functionality to test binary presence or absence rare alleles in upstream sequence
library(reshape)
library(ggplot2)
library(gridExtra)
smile_plotter=function(count_table, count_col, model_p, main, col, ylab){
  model = lm(count_table[,count_col] ~ poly(count_table$Expr_bin, 2))
  data=data.frame(expr=seq(1, dim(count_table)[1], by=1))
  data$yhat=predict(model, data) 
  R2=round(summary(model)$r.squared, digits=3)
  f=summary(model)$fstatistic
  model_p=as.character(pf(f[1],f[2],f[3], lower.tail=F))
  pval1=substring(text=model_p, first=0, last=4)
  pval2=substring(text=model_p, first=gregexpr(pattern="e", model_p), last=nchar(model_p))
  predicted.intervals <- predict(model,data.frame(x=unique(count_table$Expr_bin)),interval='confidence', level=0.999)
  upper=predicted.intervals[,"lwr"]
  lower=predicted.intervals[,"upr"]
  #return(summary(model))
  # if (max(count_table[,count_col])<=1.3){
  #   pdf(paste(basedir, "plots/", main, " std_y_ax.pdf", sep=""), height=10, width=10)
  #   plot(count_table[,count_col] ~ count_table$Expr_bin, ylim=c(0, 1.3), pch=19, cex=1.5, cex.axis=1.65, cex.lab=1.6, col="deepskyblue4", xlab="Expression rank", ylab="Fraction of ind's w/ >=1 rare variant 5kb upstream of gene", main=main)
  #   lines(x=count_table$Expr_bin, y=data$yhat, lwd=3)  
  #   polygon(c(unique(count_table$Expr_bin), rev(unique(count_table$Expr_bin))), c(lower, rev(upper)), col=adjustcolor("grey30", alpha=0.35), border=NA)
  #   legend(cex = 1.75, "top", bty="n", legend=as.expression(c(bquote(R^2 == .(R2)), paste("P = ", pval1, pval2, sep=""))))
  #   dev.off()
  # }
  #plot w/o standardized y axis
  #pdf(paste(basedir, "plots/", main, ".pdf", sep=""), height=10, width=10)
  plot(count_table[,count_col] ~ count_table$Expr_bin, pch=19, cex=1.5, cex.axis=1.65, cex.lab=1.6, col=col, xlab="Expression rank", ylab=ylab, main=main) 
  #plot(count_table[,count_col] ~ count_table$Expr_bin, pch=19, cex=1.5, cex.axis=1.65, cex.lab=1.6, col=col, xlab="Expression rank", ylab=expression("Mean ct of rare SNPs 5kb upstr. of gene"), main=main) 
  
  lines(x=count_table$Expr_bin, y=data$yhat, lwd=3) 
  polygon(c(unique(count_table$Expr_bin), rev(unique(count_table$Expr_bin))), c(lower, rev(upper)), col=adjustcolor("grey30", alpha=0.35), border=NA)
  legend(cex = 1.75, "top", bty="n", legend=as.expression(c(bquote(R^2 == .(R2)), paste("P = ", pval1, pval2, sep=""))))
  #dev.off()
}
#The following is 
all_gene_5kb_rare_ct_orig=read.table(header=T, file="/media/kak268/A_2TB_Internal/RNA/Expressions_quants/HTSEQ_counts_STAR_against_Zea_mays_B73_AGPv3.29_w_gtf_annotation/count_mat/count_mats_matched_to_AltNames/expr_vs_rare_alleles_AvgOverDups/Upstream_rare_allele_counts_with_expr_bins/all_genes/5kb_rare_ct_all_genes_w_hmp321_taxa_all_parts.txt")

top_gene_number_blocks=c(1) #OR 5001
Orig_rare_count_Or_rare_pres_binary=1 # binary refers to just counting if presence of at least one rare allele correlates with extreme expression, which was done in response to a reviewer
tissue_list=c("L3Base", "L3Tip", "LMAD", "LMAN", "Kern", "GRoot", "GShoot")
basedir="/media/kak268/A_2TB_Internal/RNA/Expressions_quants/HTSEQ_counts_STAR_against_Zea_mays_B73_AGPv3.29_w_gtf_annotation/count_mat/count_mats_matched_to_AltNames/all_non0_median_exp_genes_by_tissue_AvgOverDups/"

gene_block_size=4999
if (Orig_rare_count_Or_rare_pres_binary==1){
  all_gene_5kb_rare_ct_outside_loop=all_gene_5kb_rare_ct_orig
  colnames(all_gene_5kb_rare_ct_outside_loop)[3]="Upstr_5kb_rare_ct"
  color="turquoise3"
  cat("color is ", color)
  plotname=paste("Expression vs 5kb upstream homozy. rare allele (MAF<0.05) count \n for top", top_gene_number_blocks," to ", top_gene_number_blocks+gene_block_size ," genes with highest median exp", sep="")
  plotname_w_o_subset=paste("Expression vs 5kb upstream homozy. rare allele (MAF<0.05) count for top ",top_gene_number_blocks[1], " to ",top_gene_number_blocks[1]+4999 , " genes based on highest med exp", sep="")
  ylab="Mean ct of rare SNPs 5kb upstr. of gene"
}else if(Orig_rare_count_Or_rare_pres_binary==2){
  all_gene_5kb_rare_ct_outside_loop=all_gene_5kb_rare_ct_orig
  all_gene_5kb_rare_ct_outside_loop[all_gene_5kb_rare_ct_outside_loop$Upstr_5kb_rare_homo_ct>=1,3]=1
  colnames(all_gene_5kb_rare_ct_outside_loop)[3]="Upstr_5kb_rare_ct"
  plotname=paste("Expression vs 5kb upstream homozy. rare allele (MAF<0.05) pres. abs. \n count for top", top_gene_number_blocks," to ", top_gene_number_blocks+gene_block_size ," genes with highest median exp", sep="")
  plotname_w_o_subset=paste("Expression vs fraction with at least 1 5kb upstream homozy. rare allele (MAF<0.05) for top ",top_gene_number_blocks[1], " to ",top_gene_number_blocks[1]+4999 , " genes based on highest med exp", sep="")
  color="orange4"
  cat("color is ", color)
  ylab=paste("Fraction of ind's w/ ",expression(">=")," 1 rare SNP 5kb upstr. of gene", sep="")
}
pdf(paste(basedir, "plots/smiles/", plotname_w_o_subset, "smileplot.pdf", sep=""), height=20, width=20)
mar.default=c(5,4,4,1)
par(mfrow=c(3,3),mar = (mar.default + c(0, 1, 0, 0)))

write.table(x = "Top gene fracts divided by headers below",file=paste(basedir, "plots/smiles/", plotname_w_o_subset, ".txt", sep=""), col.names = F, row.names = F, quote = F)
for (h in tissue_list){
  cat("tissue is ", h)
  taxa_exp_counts=read.table(file=paste(basedir, h, "_all_non0_expressed_genes_w_STAR_HTSEQ_DESEQ_actual_counts.txt", sep=""), header=T)
  gene_median_rank=read.table(file=paste(basedir, h, "_all_non0_expressed_genes_w_STAR_HTSEQ_DESEQ_median_counts.txt", sep=""), header=T)
  table_to_rec_pvals=NULL
  #setup output file for data underlying the plots
  
  for(q in 1:length(top_gene_number_blocks)){
    top_number=top_gene_number_blocks[q]
    cat("lower edge of gen block is ", top_number)
    plotname=paste("Expression vs 5kb upstream homozy. rare allele (MAF<0.05) count \n for top ",top_number," to ", top_number+gene_block_size ," genes based on highest med exp ",h, sep="")
    data_top_geneset=gene_median_rank[top_number:(top_number+gene_block_size),1]
    #data_top_geneset <- data.frame(lapply(data_top_geneset, as.character), stringsAsFactors=FALSE)
    
    #keep rare cts only for specified genes
    all_gene_5kb_rare_ct=all_gene_5kb_rare_ct_outside_loop[all_gene_5kb_rare_ct_outside_loop$Feature_Name %in% data_top_geneset,]
    
    #keep expr values only or specified genes
    top_genes_w_exr_values = taxa_exp_counts[, c("HMP32Name",as.character(data_top_geneset))]
    bin_top_gene_w_exr_values=cbind(Taxon=top_genes_w_exr_values$HMP32Name, as.data.frame(apply(top_genes_w_exr_values[,2:dim(top_genes_w_exr_values)[2]], rank, MARGIN=2, ties.method="random")))
    melt_expr_bins=melt(bin_top_gene_w_exr_values,variable_name="Taxon")
    colnames(melt_expr_bins)[2]="Feature_Name"
    
    all_gene_5kb_rare_ct_top_X=all_gene_5kb_rare_ct[all_gene_5kb_rare_ct$Feature_Name %in% unique(melt_expr_bins$Feature_Name), ]
    
    #switch between all_gene_5kb_rare_ct_top_5000 AND all_gene_5kb_rare_ct_next_5000 and all_gene_5kb_rare_ct_all_median_Non0
    ex_bins_w_5kb_ct=merge(all_gene_5kb_rare_ct_top_X, melt_expr_bins, by=c("Feature_Name","Taxon"))
    head(ex_bins_w_5kb_ct)
    colnames(ex_bins_w_5kb_ct)[4]="Expr_bin"
    colnames(ex_bins_w_5kb_ct)[3]="Upstr_5kb_rare_ct"
    
    #get avg num of rare alleles per exp bin
    summed_rare_allele_counts=NULL# clear out any potential old one
    summed_rare_allele_counts_to_save=NULL
    summed_rare_allele_counts=cbind.data.frame(Expr_bin=sort(unique(ex_bins_w_5kb_ct$Expr_bin)), Upstr_5kb_rare_ct=NA)#, Upstr_5kb_rare_homohet_ct=NA, Upstr_5kb_rare_homo_ct_GERP_over_2=NA, Upstr_5kb_rare_het_ct_GERP_over_2=NA, Upstr_5kb_rare_homohet_ct_GERP_over_2=NA))
    for (k in summed_rare_allele_counts$Expr_bin){
      summed_rare_allele_counts[summed_rare_allele_counts$Expr_bin==k,"Upstr_5kb_rare_ct"]=sum(ex_bins_w_5kb_ct[ex_bins_w_5kb_ct$Expr_bin==k,"Upstr_5kb_rare_ct"])/length(unique(ex_bins_w_5kb_ct$Feature_Name))
    }
    summed_rare_allele_counts_to_save=summed_rare_allele_counts
    colnames(summed_rare_allele_counts_to_save)[1]=paste(top_number,"_to_", top_number+gene_block_size,"_top_genes_", h, "_Expr_bin",sep="")
    write.table(summed_rare_allele_counts_to_save, file=paste(basedir, "plots/smiles/", plotname_w_o_subset, ".txt", sep=""), row.names = F, quote = F, append=T)
    
    #plot smiles
    
    #underlying_data_for_plot=cbind(underlying_data_for_plot,summed_rare_allele_counts_to_save)
    smile_plotter(count_table=summed_rare_allele_counts, col=color, count_col="Upstr_5kb_rare_ct", ylab=ylab, main=paste(plotname, sep=""))
  }
  
}
dev.off()





##
##
##Plot bar plots for all tissues with different fractions of top 10k subset
library(reshape)
library(ggplot2)
library(gridExtra)
basedir="/media/kak268/A_2TB_Internal/RNA/Expressions_quants/HTSEQ_counts_STAR_against_Zea_mays_B73_AGPv3.29_w_gtf_annotation/count_mat/count_mats_matched_to_AltNames/all_non0_median_exp_genes_by_tissue_AvgOverDups/"
all_gene_5kb_rare_ct_orig=read.table(header=T, file="/media/kak268/A_2TB_Internal/RNA/Expressions_quants/HTSEQ_counts_STAR_against_Zea_mays_B73_AGPv3.29_w_gtf_annotation/count_mat/count_mats_matched_to_AltNames/expr_vs_rare_alleles_AvgOverDups/Upstream_rare_allele_counts_with_expr_bins/all_genes/5kb_rare_ct_all_genes_w_hmp321_taxa_all_parts.txt")

h="L3Base"#c("L3Base", "L3Tip", "LMAD", "LMAN", "Kern", "GRoot", "GShoot")
h="L3Tip"
h="LMAD"
h="LMAN"
h="Kern"
h="GRoot"
h="GShoot"
num_of_subsets=10
#for(num_of_subsets in c(2,5,10)){
if(num_of_subsets==2){
  top_gene_number_blocks=c(1,5001)
  gene_block_size=4999
}else if(num_of_subsets==5){
  top_gene_number_blocks=c(1,2001,4001,6001,8001)
  gene_block_size=1999
}else if(num_of_subsets==10){
  top_gene_number_blocks=c(1,1001,2001,3001,4001,5001,6001,7001,8001,9001)
  gene_block_size=999
}
plotdir=paste("plots/bar_plots/",length(top_gene_number_blocks),"_subsets/", sep="")
dir.create(file.path(basedir, plotdir), showWarnings = F)
tissue_list=c("L3Base", "L3Tip", "LMAD", "LMAN", "Kern", "GRoot", "GShoot")
#Cant use the loop, the plots come out empty; ggplot objects not being retained through the loops
#for (h in tissue_list){ 
plotname_w_o_subset=paste("Expression vs 5kb upstream homozy. rare allele (MAF<0.05) count for ",length(top_gene_number_blocks) ," top subsets of genes based on highest med exp ",h, sep="")
#mar.default=c(5,4,4,1)
#par(mfrow=c(2,5),mar = (mar.default + c(0, 1, 0, 0)))

plot_list=list()

write.table(x = "Top gene fracts divided by headers below",file=paste(basedir, plotdir, plotname_w_o_subset, ".txt", sep=""), col.names = F, row.names = F, quote = F)

cat("tissue is ", h)
taxa_exp_counts=read.table(file=paste(basedir, h, "_all_non0_expressed_genes_w_STAR_HTSEQ_DESEQ_actual_counts.txt", sep=""), header=T)
gene_median_rank=read.table(file=paste(basedir, h, "_all_non0_expressed_genes_w_STAR_HTSEQ_DESEQ_median_counts.txt", sep=""), header=T)
table_to_rec_pvals=NULL
#setup output file for data underlying the plots

for(q in 1:length(top_gene_number_blocks)){
  top_number=top_gene_number_blocks[q]
  cat("lower edge of gen block is ", top_number)
  plotname=paste("Expression vs 5kb upstream homozy. rare allele (MAF<0.05) count \n for top ",top_number," to ", top_number+gene_block_size ," genes based on highest med exp ",h, sep="")
  data_top_geneset=gene_median_rank[top_number:(top_number+gene_block_size),1]
  #data_top_geneset <- data.frame(lapply(data_top_geneset, as.character), stringsAsFactors=FALSE)
  
  #reset
  all_gene_5kb_rare_ct=all_gene_5kb_rare_ct_orig
  #keep rare cts only for specified genes
  all_gene_5kb_rare_ct=all_gene_5kb_rare_ct[all_gene_5kb_rare_ct$Feature_Name %in% data_top_geneset,]
  
  #keep expr values only or specified genes
  top_genes_w_exr_values=taxa_exp_counts[, c("HMP32Name",as.character(data_top_geneset))]
  bin_top_gene_w_exr_values=cbind(Taxon=top_genes_w_exr_values$HMP32Name, as.data.frame(apply(top_genes_w_exr_values[,2:dim(top_genes_w_exr_values)[2]], rank, MARGIN=2, ties.method="random")))
  melt_expr_bins=melt(bin_top_gene_w_exr_values,variable_name="Taxon")
  colnames(melt_expr_bins)[2]="Feature_Name"
  
  all_gene_5kb_rare_ct_top_X=all_gene_5kb_rare_ct[all_gene_5kb_rare_ct$Feature_Name %in% unique(melt_expr_bins$Feature_Name), ]
  
  #switch between all_gene_5kb_rare_ct_top_5000 AND all_gene_5kb_rare_ct_next_5000 and all_gene_5kb_rare_ct_all_median_Non0
  ex_bins_w_5kb_ct=merge(all_gene_5kb_rare_ct_top_X, melt_expr_bins, by=c("Feature_Name","Taxon"))
  head(ex_bins_w_5kb_ct)
  colnames(ex_bins_w_5kb_ct)[4]="Expr_bin"
  colnames(ex_bins_w_5kb_ct)[3]="Upstr_5kb_rare_ct"
  
  #get avg num of rare alleles per exp bin
  summed_rare_allele_counts=NULL# clear out any potential old one
  summed_rare_allele_counts_to_save=NULL
  summed_rare_allele_counts=cbind.data.frame(Expr_bin=sort(unique(ex_bins_w_5kb_ct$Expr_bin)), Upstr_5kb_rare_ct=NA)#, Upstr_5kb_rare_homohet_ct=NA, Upstr_5kb_rare_homo_ct_GERP_over_2=NA, Upstr_5kb_rare_het_ct_GERP_over_2=NA, Upstr_5kb_rare_homohet_ct_GERP_over_2=NA))
  for (k in summed_rare_allele_counts$Expr_bin){
    summed_rare_allele_counts[summed_rare_allele_counts$Expr_bin==k,"Upstr_5kb_rare_ct"]=sum(ex_bins_w_5kb_ct[ex_bins_w_5kb_ct$Expr_bin==k,"Upstr_5kb_rare_ct"])/length(unique(ex_bins_w_5kb_ct$Feature_Name))
  }
  summed_rare_allele_counts_to_save=summed_rare_allele_counts
  colnames(summed_rare_allele_counts_to_save)[1]=paste(top_number,"_to_", top_number+gene_block_size,"_top_genes_", h, "_Expr_bin",sep="")
  write.table(summed_rare_allele_counts_to_save, file=paste(basedir, plotdir, plotname_w_o_subset, ".txt", sep=""), row.names = F, quote = F, append=T)
  
  max_bin=max(ex_bins_w_5kb_ct$Expr_bin)
  min_bin=min(ex_bins_w_5kb_ct$Expr_bin)
  max_bins_count_values=(ex_bins_w_5kb_ct[which(ex_bins_w_5kb_ct$Expr_bin %in% c(max_bin, max_bin-1, max_bin-2, max_bin-3, max_bin-4)),"Upstr_5kb_rare_ct"])
  mid_bins_count_values=(ex_bins_w_5kb_ct[which(ex_bins_w_5kb_ct$Expr_bin %in% c(round(max_bin*.25,digits = 0):round(max_bin*0.75,digits = 0))),"Upstr_5kb_rare_ct"])
  min_bins_count_values=(ex_bins_w_5kb_ct[which(ex_bins_w_5kb_ct$Expr_bin %in% c(min_bin, min_bin+1, min_bin+2, min_bin+3, min_bin+4)),"Upstr_5kb_rare_ct"])
  
  
  i=plotname
  a=t.test(min_bins_count_values, max_bins_count_values)
  b=t.test(mid_bins_count_values, max_bins_count_values)
  c=t.test(mid_bins_count_values, min_bins_count_values)
  single_gene_set_table_to_rec_pvals=c(Tissue=h,Top_window=c(paste(top_number," to ", top_number+gene_block_size, sep = "")), HighVLow=a$p.value, HighVMid=b$p.value, LowVMid=c$p.value)
  if(!exists("table_to_rec_pvals")){
    table_to_rec_pvals=single_gene_set_table_to_rec_pvals
  }else{
    table_to_rec_pvals=rbind(table_to_rec_pvals, single_gene_set_table_to_rec_pvals)  
  }
  
  #if(!exists("pval_table")){
  pval_table=as.data.frame(cbind(Filename=basename(i), P.value=a[["p.value"]], Expression_bin_type="Low", Expression_bins_counts=min_bins_count_values))
  pval_table=rbind(pval_table, cbind(Filename=basename(i), P.value=a[["p.value"]], Expression_bin_type="Mid", Expression_bins_counts=mid_bins_count_values))
  pval_table=rbind(pval_table, cbind(Filename=basename(i), P.value=a[["p.value"]], Expression_bin_type="High", Expression_bins_counts=max_bins_count_values))
  # }else{
  #   pval_table=rbind(pval_table, cbind(Filename=basename(i), P.value=a[["p.value"]], Expression_bin_type="Low", Expression_bins_counts=min_bins_count_values))
  #   pval_table=rbind(pval_table, cbind(Filename=basename(i), P.value=a[["p.value"]], Expression_bin_type="Mid", Expression_bins_counts=mid_bins_count_values))
  #   pval_table=rbind(pval_table, cbind(Filename=basename(i), P.value=a[["p.value"]], Expression_bin_type="High", Expression_bins_counts=max_bins_count_values))
  # }
  pval_table=as.data.frame(pval_table)
  pval_table$Expression_bins_counts=as.numeric(levels(pval_table$Expression_bins_counts))[pval_table$Expression_bins_counts]
  pval_table$P.value=sprintf("%.2e",as.numeric(levels(pval_table$P.value))[pval_table$P.value])
  unique(pval_table$P.value)
  
  min.mean.sd.max <- function(x) {
    r <- c(  quantile(x,0.05), quantile(x,0.05), mean(x), quantile(x,0.95), quantile(x,0.95))
    names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
    r
  }
  
  myData <- aggregate(pval_table$Expression_bins_counts, by = list(HiLo = pval_table$Expression_bin_type, Gene_set = pval_table$Filename),
                      FUN = function(x) c(mean = mean(x), sd = sd(x), n = length(x)))
  
  #add the enrichement in each gene set over the middle two quartiles
  myData = cbind(myData, enrich_fraction_over_middle_2_quartiles=0, enrich_fraction_over_highest_5_ranks=0)
  myData <- do.call(data.frame, myData)
  for (q in unique(myData$Gene_set)){
    myData[myData$Gene_set==q,"enrich_fraction_over_middle_2_quartiles"]=myData[myData$Gene_set==q,"x.mean"]/myData[myData$Gene_set==q & myData$HiLo=="Mid","x.mean"]
    myData[myData$Gene_set==q,"enrich_fraction_over_highest_5_ranks"]=myData[myData$Gene_set==q,"x.mean"]/myData[myData$Gene_set==q & myData$HiLo=="High","x.mean"]
  }
  #overall_enrichment_mat=rbind(overall_enrichment_mat, myData)
  
  # myData$mean=myData$x[,"mean"]
  # myData$sd=myData$x[,"sd"]
  # myData$n=myData$x[,"n"]
  
  myData$se=myData$x.sd/sqrt(myData$x.n)
  myData$y_min=myData$x.mean - myData$se
  myData$y_max=myData$x.mean + myData$se
  
  #This reorders the rows so that they reflect the intended order in the plot
  myData$HiLo=factor(myData$HiLo, levels = c("Low", "Mid", "High"))
  #myData$Gene_set=factor(myData$Gene_set, levels = c(paste(h, "_0kb_to_5kbupstr_top_5000_expressed_genes_w_STAR_HTSEQ_DESEQ_rare_allele_cts.txt", sep=""), paste(h,  "_0kb_to_5kbupstr_top5001to10000_non0_expressed_genes_w_STAR_HTSEQ_DESEQ_rare_allele_cts.txt", sep=""), paste(h, "_0kb_to_5kbupstr_rand_5000_expressed_genes_w_STAR_HTSEQ_DESEQ_rare_allele_cts.txt", sep="")))
  
  myData=myData[order(myData$Gene_set, myData$HiLo),]
  write.table(myData, file=paste(basedir, plotdir, plotname, "enrichment_by_group.txt", sep=""), row.names = F, quote = F)
  
  p=NULL
  p =ggplot(data = myData, aes(x = factor(Gene_set), y = x.mean, fill = factor(HiLo))) +
    geom_bar(stat = "identity", position = position_dodge(0.9)) +
    coord_cartesian(ylim = c(0, 3.5))+ #TEST
    #Note with the current version of ggplot (aug 7 2017) the error bars for high and low bins were reversed
    geom_errorbar(ymin=rev(myData$y_min), ymax=rev(myData$y_max), position = position_dodge(0.9), width = 0.25) +
    #ggtitle(paste(h, " smile plot extreme low vs \n high bins means and se", sep="")) +
    ggtitle(plotname) +
    labs(x= "Gene set", y = paste(h," upstream rare allele count",sep = "")) +
    #scale_x_discrete(labels=c( "BUSCO conserved \nplantae orthologs",  "Random expr. genes", "Top expr. genes", "Medium expr. genes")) +
    #scale_x_discrete(labels=c("Top expr. genes", "Medium expr. genes", "Random non 0 \n expr. genes")) +
    scale_x_discrete(labels=c("")) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))+
    #theme(axis.text=element_text(size=18))
    scale_fill_manual(values=c(Low="deeppink3", High="steelblue3", Mid="goldenrod2"),guide=guide_legend(title="5 extreme expr. bins \n vs mid. two quartiles")) +
    #annotate("text",x=1,y=max(myData$y_max[1:3])*1.025,label=paste("P = ",unique(pval_table$P.value)[1], sep = ""))+
    #annotate("text",x=1,y=max(myData$y_max[1:3])*1.08,label=paste("P = ",unique(pval_table$P.value)[1], sep = ""))+
    #annotate("text",x=2,y=max(myData$y_max[4:6])*1.08,label=paste("P = ",unique(pval_table$P.value)[2], sep = ""))+
    #annotate("text",x=3,y=max(myData$y_max[7:9])*1.08,label=paste("P = ",unique(pval_table$P.value)[3], sep = ""))+
    #first_y=max(myData$y_max[1:3])*1.035
    
    #annotate("text",x=4,y=max(myData$y_max[10:12])*1.025,label=paste("P = ",unique(pval_table$P.value)[4], sep = ""))+
    theme_classic(base_size=15)#, envir = .GlobalEnv)
  #print(p)
  plot_list[[q]]=p
  
}
write.table(x=table_to_rec_pvals, file =paste(basedir,plotdir,plotname_w_o_subset, " pval_table.txt", sep=""), quote=F, row.names = F)

cat("in bar plot print function")
#https://stackoverflow.com/questions/33258000/plot-multiple-ggplot2-on-same-page
pdf(paste(basedir, plotdir, plotname_w_o_subset, "barplot.pdf", sep=""), height=6, width=6*length(top_gene_number_blocks), onefile = F)
marrangeGrob(plot_list, nrow = 1, ncol = length(top_gene_number_blocks))
#grid.arrange(plot_list, nrow = 1, ncol = length(top_gene_number_blocks),newpage = T)
dev.off()
#}

