#!/usr/bin/env Rscript

rm(list=ls())

setwd("/bioinfo1/projects/bean/VEF/genomic_selection/scripts")

args = commandArgs(trailingOnly=TRUE)
vargs <- strsplit(args, ",")

if (length(args) < 4) {
  
  stop("At least four arguments must be supplied:\n
  \n  Usage: Rscript 3_run_trait <prior1>,...,<priorN> <trait1>,...,<traitN> <phenotype(s)> <outDir>\n
  Positional arguments:
  \tprior     :  Comma separated list of priors to be tested.  Available options are 
  \t\t    BayesA,BayesB,BayesC,BayesRR,BLasso,FIXED,BLassof,RKHS,GBLUP\n 
  \ttrait     :  Comma separated list of traits to be tested\n
  \tphenotypes:  File with the phenotype values for every line. One column per trait. If two phenotypes are provided, then the first one will be used for the Training population, and the second will be used as the Breeding population.
  \t\t    First row is header line.\n
  \toutDir    :  A directory path where the output wil be stored\n
  \tcombinat  :\ This is an optional argument. In case you already have a matrix with population partitions, you can provide it in here.\n\n", call.=FALSE)
  
} else {
  
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
  if (exists('all_phen_gen2')){
    cat('\n\n',length(all_phen_gen2)-nrow(combinat),' phenotypes will be predicted from zero.\n\n')
  } 
  
  cat('model\ttrait\trandomPop\tcorr\tstartedAt\n')
  
  for (i in 1:100){
    
    for (prior in model){
      
      for (trait in traits){
        
        set.seed(1234)
        
        if (is.na(vargs[[3]][2])){
          
          if (trait %in% names(phen)){
            myCorr = round( runBGLR(y = phen, trait = trait, X = geno , pop.split = combinat[,i], model = prior )$cor, 5)
            cat(prior,"\t",trait,"\tpop",i,"\t",myCorr,"\t",paste(Sys.time(),'\n'), sep = '')
          }
          
        } else {
          
          shared_traits = intersect(names(phen), names(phen2))
          
          if (trait %in% shared_traits){
            myCorr = round( runBGLR(y = phen, trait = trait, X = geno , pop.split = combinat[,i], yBP = phen2, model = prior )$cor, 5)
            cat(prior,"\t",trait,"\tpop",i,"\t",myCorr,"\t",paste(Sys.time(),'\n'), sep = '')
          }
        }
        }
      }
    }
  }

