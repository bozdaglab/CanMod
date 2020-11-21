install.packages("bnlearn")
library(data.table)
library(pbapply)
library(doParallel)
library(parallel)
library(HDCI)
library(glmnet)
library(bnlearn)

# get regression df for each mRNA with cna, methyl, TF and miRNA included as candidate regulators
getDF = function(regression.data, mRNA, include.lncRNA = F){
  miRNA.target.interactions = regression.data$miRNA.target.interactions #TCGA, starbase & targetScan common miRNA
  miRNA.candidates = as.character(miRNA.target.interactions[miRNA.target.interactions$target == mRNA,c("miRNA")])
  miRNA.candidates = miRNA.candidates[miRNA.candidates %in% rownames(regression.data$miRNA)]
  
  tf.target.interactions = regression.data$tf.target.interactions
  TF.candidates = tf.target.interactions$tf[which(tf.target.interactions$target == mRNA)] 
  TF.candidates = setdiff(TF.candidates[TF.candidates %in% rownames(regression.data$mRNA)], mRNA)
  
  mrna.expression.vector = regression.data$mRNA[mRNA,]
  mrna.cna.vector = regression.data$cna[mRNA,]
  mrna.methyl.vector = regression.data$methyl[mRNA,]
  predicted.TF.df = regression.data$mRNA[TF.candidates,]
  predicted.miRNA.df = regression.data$miRNA[miRNA.candidates,]
  
  
  if (include.lncRNA){
    if (nrow(predicted.miRNA.df) < 1){
      return(NULL)
    }else{
      regression.df = as.data.frame(scale(cbind(t(mrna.expression.vector), t(predicted.miRNA.df))))
      colnames(regression.df)[1] = c("mRNA")
      return(regression.df)
    }
  }
  
  if (nrow(predicted.TF.df) >= 1 & nrow(predicted.miRNA.df) < 1){
    regression.df = as.data.frame(scale(cbind(t(mrna.expression.vector),t(mrna.cna.vector),t(mrna.methyl.vector),
                                              t(predicted.TF.df))))
  }else if (nrow(predicted.TF.df) < 1 & nrow(predicted.miRNA.df) >= 1){
    regression.df = as.data.frame(scale(cbind(t(mrna.expression.vector),t(mrna.cna.vector),t(mrna.methyl.vector),
                                              t(predicted.miRNA.df))))
  }else if (nrow(predicted.TF.df) >= 1 & nrow(predicted.miRNA.df) >= 1){
    regression.df = as.data.frame(scale(cbind(t(mrna.expression.vector),t(mrna.cna.vector),t(mrna.methyl.vector),
                                              t(predicted.miRNA.df),t(predicted.TF.df))))
  }else{
    regression.df = NULL
  }
  
  if (is.null(regression.df)){
    return(NULL)
  }
  
  colnames(regression.df)[1:3] = c("mRNA","CNA","Methyl")
  regression.df$CNA[which(is.na(regression.df$CNA))] = 0
  regression.df$Methyl[which(is.na(regression.df$Methyl))] = 0
  
  return(regression.df)
}

