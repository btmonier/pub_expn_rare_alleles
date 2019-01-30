#Kal Kremling
#

library(plyr)
library("foreach")
library("doMC")
#library("parallel")
registerDoMC(cores=12)
library(reshape2)
library(readr)
#install.packages("heR.misc") # for geometric mean calculations
#library(heR.misc)

#install.packages("glmnet")
library(glmnet)
#Calculate the expression RESIDUALS after factoring out Hidden Factors calculated from expression

#install.packages("caret")
library(caret)

library(stringr)
library(ggplot2)

tissue_list=c("Kern", "GRoot", "L3Base", "L3Tip", "LMAD", "LMAN", "GShoot")
#tissue_list=c("Kern", "GRoot", "L3Base")


#PHENOS
#HOLLAND 06 07 Seed Weight BLUPs
if(exists("phenos")){rm(phenos)}
#Load in phenotypes 
#phenos=read.table("/home/kak268/work_data/maize_phenotypes/Ames_BLUPs/Merged_Ames_traits", header=T)
phenos=read.table(file="/home/kak268/work_data/maize_phenotypes/Holland_NAM_and_282_BLUPs/282 BLUPs 0607 all traits.txt", sep="\t", header=T)


#Top 5000 genes from each tissue avg over dups
basedir="/media/kak268/A_2TB_Internal/RNA/Expressions_quants/HTSEQ_counts_STAR_against_Zea_mays_B73_AGPv3.29_w_gtf_annotation/count_mat/count_mats_matched_to_AltNames/Top_5000_by_tissue_AvgOverDups/"
filelist=list.files(basedir, pattern="actual_counts.txt")
#expr in original scale dir
#plot_out_dir="/home/kak268/work_data/maize_phenotypes/Yield_BLUPs_on_Testers/plots/expression_dev_v_top_5000_w_sweet_and_pop_exclusions_lasso_and_ridge_AvgOverDups/"
#scale expr 0 1 directory
plot_out_dir="/home/kak268/work_data/maize_phenotypes/Yield_BLUPs_on_Testers/plots/expression_dev_v_top_5000_w_sweet_and_pop_exclusions_lasso_and_ridge_AvgOverDups_exp_nested_CV_10reps_no_popsweettrop/"
#plot_out_dir="/home/kak268/work_data/maize_phenotypes/Yield_BLUPs_on_Testers/plots/expression_dev_v_top_5000_w_sweet_and_pop_exclusions_lasso_and_ridge_AvgOverDups_exp_norm_to_0_1_nested_CV_10reps_no_popsweettrop/"
exclude=T
#exclusion_list=c("282set_Ia5125", "282set_Il101", "282set_Il14H", "282set_P39Goodman-Buckler", "282set_i1677a", "282set_IA2132Goodman-Buckler")
#pop and sweet exc list
exclusion_list=c("282set_Ia5125", "282set_Il101", "282set_Il14H", "282set_P39Goodman-Buckler", "282set_i1677a", "282set_IA2132Goodman-Buckler", "282set_IDS91", "282set_4722", "282set_HP301", "282set_I29", "282set_IDS28", "282set_IDS69", "282set_IDS91", "282set_SA24", "282set_Sg1533", "282set_Sg18")
filesuffix="in 5k most exp genes 0notSetToNA_AvgOverDups no sweet no pop glmmnet ridge "
use_residuals=F
avg_over_dups=T
avg_CV_for_plot=T
drop_tropicals=T
scramble_phenos=F
scale_0_1=T
#lambdas for cv.glmnet, no longer hard coding
#lambda_range=((1:40) * 0.05)

#cross validation params for Caret:
#eGrid <- expand.grid(.alpha = 0,  .lambda = (1:40) * 0.05)
#set up the 20 reps of 5-fold CV with caret https://machinelearningmastery.com/how-to-estimate-model-accuracy-in-r-using-the-caret-package/
#model predictions are kepth for each of the 20 repetitions of the 5 fold (80-20) split.
#Control <- trainControl(method = "repeatedcv", number=5, repeats = 20, verboseIter =TRUE, allowParallel = T, savePredictions = "final" )

#Reps of outer cross validation
repititions=10
#number of outer and inner cross validation folds
K_outer <- 10
K_inner <- 10

filelist=filelist[-5]#this is to remove the "L3Mid" tissue
for (i in filelist){
      cat(i, "\n")
      print(Sys.time())
      if (i=="L3Mid_top_5000_expressed_genes_w_STAR_HTSEQ_DESEQ_actual_counts.txt") {
        next
      }
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
      }
      
      if(scale_0_1==T){
        gene_w_exr_values_0_1=as.data.frame(apply(as.matrix((gene_w_exr_values[,2:dim(gene_w_exr_values)[2]])), MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X))))
        gene_w_exr_values_0_1=cbind(gene_w_exr_values$HMP32Name, gene_w_exr_values_0_1)
        gene_w_exr_values=gene_w_exr_values_0_1
        colnames(gene_w_exr_values)[1]="HMP32Name"
      }
      
      
      means=colMeans(gene_w_exr_values[,2:ncol(gene_w_exr_values)])
      var_mat=gene_w_exr_values
      var_mat[,-1]=NA # make an empty mat of same dimensions, leave names in first column
      for (j in 1:nrow(gene_w_exr_values)){
        var_mat[j,2:ncol(var_mat)]=gene_w_exr_values[j,2:ncol(var_mat)]-means
        cat(j)
      }
      var_mat[,2:ncol(var_mat)]=abs(var_mat[,2:ncol(var_mat)])
      #Load in name converter so names can be converted to Ames names
      name_decoder=read.table(header=T,"/media/kak268/A_2TB_Internal/RNA/Expressions_quants/renaming_scratch_folder/orig_1960names_and_alternate_names_w_subpops.txt")  
      name_decoder_tissue=name_decoder[name_decoder$TissueWODate==tissue_name,]
      ###Make the decoder unique
      unique_name_decoder_tissue=name_decoder_tissue[!duplicated(name_decoder_tissue$HMP32Name),]
      unique_name_decoder_tissue=unique_name_decoder_tissue[,c("HMP32Name", "RawPhenotypeNames", "Subpopulation")]
      #remove mean_dev_by_taxon_w_names=merge(gene_w_exr_values, unique_name_decoder_tissue, by="HMP32Name")
    
    #PHENOS
    if(exists("phenos_w_alt_names")){phenos_w_alt_names=NULL}
    phenos_w_alt_names=merge(unique_name_decoder_tissue, phenos, by.x="RawPhenotypeNames", by.y="TaxaName")
    if(drop_tropicals==T){
      phenos_w_alt_names=phenos_w_alt_names[!phenos_w_alt_names$Subpopulation %in% c("ts", "landrace_or_CIMMYT"),]
    }
    seed_blups_w_hmp3_name=phenos_w_alt_names[,c("HMP32Name", "TotalKernelWeight0607summer_blup")]
    seed_blups_w_hmp3_name=seed_blups_w_hmp3_name[!is.na(seed_blups_w_hmp3_name$TotalKernelWeight0607summer_blup),]
    #note the line below is a permutation test with scrambled phenotypes
    if(scramble_phenos==T){
    seed_blups_w_hmp3_name$TotalKernelWeight0607summer_blup=sample(seed_blups_w_hmp3_name$TotalKernelWeight0607summer_blup, replace = F)
    }
    
    
    
    ##Predict and plot on EXPR VALUES
    pdf(file = paste(plot_out_dir,filesuffix," ",tissue_name, " plot each cross valid vs_0607_seed_wt_BLUP_fulltable.pdf", sep=""), height = 9.2, width = 16)
    par(mfrow=c(1,2))
    
    if(exists("expr_vals_and_blups")){rm(expr_vals_and_blups)}
    expr_vals_and_blups=merge(x=seed_blups_w_hmp3_name, gene_w_exr_values,  by="HMP32Name")
    x=expr_vals_and_blups[,3:dim(expr_vals_and_blups)[2]]
    y=expr_vals_and_blups[,"TotalKernelWeight0607summer_blup"]
    
    #netFit_exp_vals <- train(x =as.matrix(x), y = y,
    #                         method = "glmnet",
    #                         tuneGrid = eGrid,
    #                         trControl = Control)
    
    
    y_obs_and_pred_expr_vals <- NULL
    for(rep_num in 1:repititions){
      #Outer and inner cross folds
      
      #Make a vector of cross fold membership to be used in subsetting data
      #fix to be flexible with non-multiples of K
      #take only the first 1:length of y to make sure you don't have more samples in the fold than actual samples
      fold <- sample(rep(x = 1:K_outer, each=ceiling(length(y)/K_outer)))[1:length(y)]
      
      #outer cross validation
      for (k in 1:K_outer) {
        y_train <- y
        y_train[fold == k] <- NA
        
        len_temp_df=length(y[fold == k])
        temp_y_obs_and_pred_expr_vals <- data.frame(y_obs=double(len_temp_df),
                                          y_pred=double(len_temp_df),
                                          Rep=integer(len_temp_df),
                                          Fold=integer(len_temp_df),
                                          lambda.min=double(len_temp_df),
                                          rowIndex=integer(len_temp_df),
                                          stringsAsFactors=FALSE)
        
        #inner cross validation w/in cv.glmnet for choosing lambda and betas Ridge Regression
        CV <- cv.glmnet(x = as.matrix(x)[!is.na(y_train), ], y = y_train[!is.na(y_train)], alpha=0, nfolds = K_inner)
        
        #y_pred[fold == k] <- predict(CV, newx = as.matrix(x)[fold == k, ], s=c("lambda.min"))[,1]
        
        temp_y_obs_and_pred_expr_vals$y_obs=y[fold == k]
        #fit the model to unseen x
        temp_y_obs_and_pred_expr_vals$y_pred=predict(CV, newx = as.matrix(x)[fold == k,], s=c("lambda.min"))[,1]
        temp_y_obs_and_pred_expr_vals$lambda.min=CV$lambda.min
        temp_y_obs_and_pred_expr_vals$Rep=rep_num
        temp_y_obs_and_pred_expr_vals$Fold=k
        temp_y_obs_and_pred_expr_vals$rowIndex=which(fold %in% k)
        
        
        if(!exists("y_obs_and_pred_expr_vals")){
          y_obs_and_pred_expr_vals=temp_y_obs_and_pred_expr_vals
        }else{
          y_obs_and_pred_expr_vals=rbind(y_obs_and_pred_expr_vals, temp_y_obs_and_pred_expr_vals)
        }
      }
    }
    
    
    plot(main=paste(tissue_name, "pred on top 5k expr genes", sep=""), y_obs_and_pred_expr_vals$y_obs, y_obs_and_pred_expr_vals$y_pred, xlab="Seed Wt. BLUPs", ylab="Mean of Pred. Seed Wt.")
    abline(fit_vals <- lm(y_obs_and_pred_expr_vals$y_pred ~ y_obs_and_pred_expr_vals$y_obs), col='red')
    #legend("topleft", bty="n", legend=paste("r = ", format(sqrt(summary(fit_vals)$r.squared), digits=4), "\nP=",format(summary(fit_vals)$coefficients[,4][2], digits=4)))
    rm(fit_vals)
    
    #CONTINUE HERE
    #sep_fold_and_rep_vals=str_split_fixed(netFit_exp_vals$pred$Resample, "\\.", 2)
    #colnames(sep_fold_and_rep_vals)=c("Fold","Rep")
    #sep_fold_and_rep_vals_w_cor=cbind(netFit_exp_vals$pred, sep_fold_and_rep_vals)
    y_obs_and_pred_expr_vals$correl=NA
    for(local_rep in unique(y_obs_and_pred_expr_vals$Rep)){
      y_obs_and_pred_expr_vals[y_obs_and_pred_expr_vals$Rep==local_rep,]$correl=cor(y_obs_and_pred_expr_vals[y_obs_and_pred_expr_vals$Rep==local_rep,]$y_pred,y_obs_and_pred_expr_vals[y_obs_and_pred_expr_vals$Rep==local_rep,]$y_obs)
    }
    write.table(y_obs_and_pred_expr_vals , quote=F, row.names = F, col.names = T, file = paste(plot_out_dir,filesuffix," ",tissue_name, " expr vals vs_0607_seed_wt_BLUP_fulltable.txt", sep=""))

    ##Predict and plot on EXPR DEVIATIONS
    #if(exists("expr_devs_and_blups")){rm(expr_devs_and_blups)}
    expr_devs_and_blups=merge(x=seed_blups_w_hmp3_name, var_mat,  by="HMP32Name")
    
    x=NULL
    y=NULL
    x=expr_devs_and_blups[,3:dim(expr_devs_and_blups)[2]]
    y=expr_devs_and_blups[,"TotalKernelWeight0607summer_blup"]
    
    #netFit_exp_devs <- train(x =as.matrix(x), y = y,
    #                         method = "glmnet",
    #                         tuneGrid = eGrid,
    #                         trControl = Control)
    
    

    y_obs_and_pred_expr_devs <- NULL
    for(rep_num in 1:repititions){
      #Outer and inner cross folds
      
      #Make a vector of cross fold membership to be used in subsetting data
      #fix to be flexible with non-multiples of K
      #take only the first 1:length of y to make sure you don't have more samples in the fold than actual samples
      fold <- sample(rep(x = 1:K_outer, each=ceiling(length(y)/K_outer)))[1:length(y)]
      
      #outer cross validation
      for (k in 1:K_outer) {
        y_train <- y
        y_train[fold == k] <- NA
        
        len_temp_df=length(y[fold == k])
        temp_y_obs_and_pred_expr_devs <- data.frame(y_obs=double(len_temp_df),
                                                    y_pred=double(len_temp_df),
                                                    Rep=integer(len_temp_df),
                                                    Fold=integer(len_temp_df),
                                                    lambda.min=double(len_temp_df),
                                                    rowIndex=integer(len_temp_df),
                                                    stringsAsFactors=FALSE)
        
        #inner cross validation w/in cv.glmnet for choosing lambda and betas Ridge Regression
        CV <- cv.glmnet(x = as.matrix(x)[!is.na(y_train), ], y = y_train[!is.na(y_train)], alpha=0, nfolds = K_inner)
        
        #y_pred[fold == k] <- predict(CV, newx = as.matrix(x)[fold == k, ], s=c("lambda.min"))[,1]
        
        temp_y_obs_and_pred_expr_devs$y_obs=y[fold == k]
        #fit the model to unseen x
        temp_y_obs_and_pred_expr_devs$y_pred=predict(CV, newx = as.matrix(x)[fold == k,], s=c("lambda.min"))[,1]
        temp_y_obs_and_pred_expr_devs$lambda.min=CV$lambda.min
        temp_y_obs_and_pred_expr_devs$Rep=rep_num
        temp_y_obs_and_pred_expr_devs$Fold=k
        temp_y_obs_and_pred_expr_devs$rowIndex=which(fold %in% k)
        
        if(!exists("y_obs_and_pred_expr_devs")){
          y_obs_and_pred_expr_devs=temp_y_obs_and_pred_expr_devs
        }else{
          y_obs_and_pred_expr_devs=rbind(y_obs_and_pred_expr_devs, temp_y_obs_and_pred_expr_devs)
        }
      }
    }
    
    
    plot(main=paste(tissue_name, "pred on top 5k exp genes, abs. deviations from the mean", sep=""), y_obs_and_pred_expr_devs$y_obs, y_obs_and_pred_expr_devs$y_pred, xlab="Seed Wt. BLUPs", ylab="Mean of Pred. Seed Wt.")
    abline(fit_vals <- lm(y_obs_and_pred_expr_devs$y_pred ~ y_obs_and_pred_expr_devs$y_obs), col='red')
    #legend("topleft", bty="n", legend=paste("r = ", format(sqrt(summary(fit_vals)$r.squared), digits=4), "\nP=",format(summary(fit_vals)$coefficients[,4][2], digits=4)))
    rm(fit_vals)
    
    #insert here
    y_obs_and_pred_expr_devs$correl=NA
    for(local_rep in unique(y_obs_and_pred_expr_devs$Rep)){
      y_obs_and_pred_expr_devs[y_obs_and_pred_expr_devs$Rep==local_rep,]$correl=cor(y_obs_and_pred_expr_devs[y_obs_and_pred_expr_devs$Rep==local_rep,]$y_pred,y_obs_and_pred_expr_devs[y_obs_and_pred_expr_devs$Rep==local_rep,]$y_obs)
    }
    write.table(y_obs_and_pred_expr_devs , quote=F, row.names = F, col.names = T, file = paste(plot_out_dir,filesuffix," ",tissue_name, " expr devs vs_0607_seed_wt_BLUP_fulltable.txt", sep=""))
    
    
    dev.off()
    
      if(avg_CV_for_plot==T){

        pdf(file = paste(plot_out_dir,filesuffix," ",tissue_name, " plot avg. of cross valid vs_0607_seed_wt_BLUP_fulltable.pdf", sep=""), height = 9.2, width = 16)
        par(mfrow=c(1,2))


        netFit_exp_vals_avg=ddply(y_obs_and_pred_expr_vals[,c("y_pred", "y_obs", "rowIndex")], "rowIndex", numcolwise(mean), .parallel = T)
        
        plot(pch=19, main=paste(tissue_name, " pred on top 5k expr vals", sep=""), netFit_exp_vals_avg[, "y_obs"], netFit_exp_vals_avg[, "y_pred"],xlab="Seed Wt. BLUPs", ylab="Mean of Pred. Seed Wt.")
        abline(fit_vals <- lm(netFit_exp_vals_avg[, "y_pred"] ~ netFit_exp_vals_avg[, "y_obs"]), col='red')
        legend("topleft", cex=1.5, bty="n", legend=paste("r = ", format(sqrt(summary(fit_vals)$r.squared), digits=4), "\nP=",format(summary(fit_vals)$coefficients[,4][2], digits=4)))
        rm(fit_vals)

        #assign(paste(tissue_name, "_expr_vals_caretFit", sep=""), netFit_exp_vals_avg$pred)


        netFit_exp_devs_avg=ddply(y_obs_and_pred_expr_devs[,c("y_pred", "y_obs", "rowIndex")], "rowIndex", numcolwise(mean), .parallel = T)

        plot(pch=19, main=paste(tissue_name, " pred on top 5k expr gene abs. dev's from mean", sep=""), netFit_exp_devs_avg[, "y_obs"], netFit_exp_devs_avg[, "y_pred"], xlab="Seed Wt. BLUPs", ylab="Mean of Pred. Seed Wt.")
        abline(fit_devs <- lm(netFit_exp_devs_avg[, "y_pred"] ~ netFit_exp_devs_avg[, "y_obs"]), col='red')
        legend("topleft", cex=1.5, bty="n", legend=paste("r = ", format(sqrt(summary(fit_devs)$r.squared), digits=4), "\nP=",format(summary(fit_devs)$coefficients[,4][2], digits=4)))
        rm(fit_devs)

        dev.off()
      }
}



#Loop to calculate ranges of prediction accuracy
library(RColorBrewer)


#acuracy_master_table builder, can be for expression of deviation of expression depending on the directory it's pointed at.
devs_or_vals="devs"
accuracy_master_table=NULL
#table_dir="/home/kak268/work_data/maize_phenotypes/Yield_BLUPs_on_Testers/plots/expression_dev_v_top_5000_w_sweet_and_pop_exclusions_lasso_and_ridge_AvgOverDups_exp_norm_to_0_1_nested_CV/"
table_dir=plot_out_dir
for (tissue in c("GRoot", "GShoot","L3Base", "L3Tip", "LMAD", "LMAN", "Kern")){
  #ecl sweet and pop
  table=read.table(stringsAsFactors = F,header=T, paste(table_dir, "in 5k most exp genes 0notSetToNA_AvgOverDups no sweet no pop glmmnet ridge  ",tissue," expr ",devs_or_vals," vs_0607_seed_wt_BLUP_fulltable.txt", sep=""))
  #also exclude tropicals
  #table=read.table(stringsAsFactors = F,header=T, paste(plot_out_dir, "in 5k most exp genes 0notSetToNA_AvgOverDups no sweet no pop no trop glmmnet ridge  ",tissue," expr ",devs_or_vals," vs_0607_seed_wt_BLUP_fulltable.txt", sep=""))
  table=cbind.data.frame(Tissue_name=tissue, Expr_value_or_deviation=devs_or_vals ,Prediction_r = table$correl)
  table=unique(table)
  cat(tissue, "\n")
  print(range(table$Prediction_r))
  if(!exists("accuracy_master_table")){
    accuracy_master_table=table
  }else{
    accuracy_master_table=rbind.data.frame(accuracy_master_table, table)
  }
}

pdf(width = 13.29, height = 4.82, file = paste0(table_dir, "/ridge regression 10 reps 10 fols outer 10 folds inner preds from devs.pdf"))
p=ggplot(accuracy_master_table, aes(Tissue_name,  Prediction_r))
p + geom_boxplot(aes(fill=factor(Tissue_name)))+ 
  theme_classic() +
  ggtitle(paste("Ridge Regression fitness predictions using expression ", devs_or_vals, sep="")) +
  #scale_y_continuous(limits = c(0, 0.52))+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"))+
  scale_fill_brewer(palette = "RdYlBu")
dev.off()


#plot both devs and vals together after including them in the same accuracy_master_table
#to do this run the acuracy_master_table builder again with devs_or_vals set to the other option
pdf(width = 13.29, height = 4.82, file = paste0(table_dir, "/ridge regression 10 reps 10 fols outer 10 folds inner preds from vals and devs.pdf"))
devs_or_vals="vals"
#table_dir="/home/kak268/work_data/maize_phenotypes/Yield_BLUPs_on_Testers/plots/expression_dev_v_top_5000_w_sweet_and_pop_exclusions_lasso_and_ridge_AvgOverDups_exp_norm_to_0_1_nested_CV/"
table_dir=plot_out_dir
for (tissue in c("GRoot", "GShoot","L3Base", "L3Tip", "LMAD", "LMAN", "Kern")){
  #ecl sweet and pop
  table=read.table(stringsAsFactors = F,header=T, paste(table_dir, "in 5k most exp genes 0notSetToNA_AvgOverDups no sweet no pop glmmnet ridge  ",tissue," expr ",devs_or_vals," vs_0607_seed_wt_BLUP_fulltable.txt", sep=""))
  #also exclude tropicals
  #table=read.table(stringsAsFactors = F,header=T, paste(plot_out_dir, "in 5k most exp genes 0notSetToNA_AvgOverDups no sweet no pop no trop glmmnet ridge  ",tissue," expr ",devs_or_vals," vs_0607_seed_wt_BLUP_fulltable.txt", sep=""))
  table=cbind.data.frame(Tissue_name=tissue, Expr_value_or_deviation=devs_or_vals ,Prediction_r = table$correl)
  table=unique(table)
  cat(tissue, "\n")
  print(range(table$Prediction_r))
  if(!exists("accuracy_master_table")){
    accuracy_master_table=table
  }else{
    accuracy_master_table=rbind.data.frame(accuracy_master_table, table)
  }
}
p=ggplot(accuracy_master_table, aes(Tissue_name,  Prediction_r))
p + geom_boxplot(aes(fill=factor(Expr_value_or_deviation)))+
  theme_classic() +
    ggtitle(paste("Ridge Regression fitness predictions using expression devs and vals", sep=""))+ 
    #scale_y_continuous(limits = c(0, 0.52))+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"))+
    scale_fill_brewer(palette = "RdYlBu")
dev.off()

