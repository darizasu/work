#!/usr/bin/env Rscript

rm(list=ls())

setwd("/bioinfo1/projects/bean/VEF/genomic_selection/scripts")

args = commandArgs(trailingOnly=TRUE)
vargs <- strsplit(args, ",")

if (length(args) < 4) {
  
  stop("At least four arguments must be supplied:\n
  \n  Usage: Rscript 3_run_trait <prior1>,...,<priorN> <trait1>,...,<traitN> <phenotypes> <outDir>\n
  Positional arguments:
  \tprior     :  Comma separated list of priors to be tested.  Available options are 
  \t\t    BayesA,BayesB,BayesC,BayesRR,BLasso,FIXED,BLassof,RKHS,GBLUP\n 
  \ttrait     :  Comma separated list of traits to be tested\n
  \tphenotypes:  File with the phenotype values for every line. One column per trait. If two phenotypes are provided, then the first one will be used for the Training population, and the second will be used as the Breeding population. In such case, no iterations will be performed.
  \t\t    First row is header line.\n
  \toutDir   :   A directory path where the output wil be stored\n\n", call.=FALSE)
  
} else {
  
  # load models
  library(BGLR)
  
  model  = vargs[[1]]
  traits = vargs[[2]]
  phen   = vargs[[3]][1]
  phen2  = vargs[[3]][2]
  outDir = vargs[[4]]
  combinat = if (length(args) > 4) vargs[[5]] else NA
  
  source("./1_getting_data.R")
  source("./2_prepare_models.R")
  # model= c( "BayesA", "BayesB", "BayesC", "BLasso", "BLassof", "RKHS", "GBLUP" )
  
  cat('\n\n',length(all_phen_gen),'individuals will be used in this run.\n\n')
  
  if (is.na(vargs[[3]][2])){
    cat('model\ttrait\trandomPop\tstartedAt\n')
    for (i in 1:100){
      for (prior in model){
        for (trait in traits){
          cat(prior,"\t",trait,"\tpop",i,"\t",paste(Sys.time(),'\n'), sep = '')
          set.seed(1234)
          cors[[trait]][i,prior] = round( runBGLR(y = phen, trait = trait, X = geno , pop.split = combinat[,i], model = prior )$cor, 5)
          write.table(cors[[trait]], paste(outDir,"/cors_",trait,".txt",sep=''), sep="\t", row.names=FALSE, quote = F)
        }
      }
    }
  } else {
    
    shared_traits = intersect(names(phen), names(phen2))
    traits = traits[traits %in% shared_traits]
    cat('model\ttrait\tbreedingPop\ttrainingPop\tstartedAt\n')
    for(trait in traits){
      for (prior in model){
        
        
        cat(prior,"\t",trait,"\t",vargs[[3]][2],"\t",vargs[[3]][1],"\t",paste(Sys.time(),'\n'), sep = '')
        cors[[trait]][,prior] = round( runBGLR(y = phen, trait = trait, X = geno , yBP = phen2, model = prior )$cor, 5)
        write.table(cors[[trait]], paste(outDir,"/cors_",trait,".txt",sep=''), sep="\t", row.names=FALSE, quote = F)
        
      }
    }
    
  }
  
}

