#!/home/dariza/bin/Rscript

rm(list=ls())

list.of.packages = c("argparse","BGLR")
new.pckgs = list.of.packages[!(list.of.packages %in% rownames(installed.packages()))]

if (length(new.pckgs)){
  cat("\nWarning message:\nThis script will try to install the packages",new.pckgs,"and all its dependencies. Otherwise abort this execution\n")
  Sys.sleep(10)
  install.packages(pkgs=new.pckgs, repos="http://cran.r-project.org", dependencies = T)
}

library(argparse)

parser = ArgumentParser()

parser$add_argument("-m", type='character', metavar='Model(s)',
                    help="Comma separated list of priors to be tested. Available options are BayesA,BayesB,BayesC,BayesRR,BLasso,FIXED,BLassof,RKHS,GBLUP")
parser$add_argument("-t", type='character', metavar='Trait(s)',
                    help="Comma separated list of traits to be tested")
parser$add_argument("-p", type="character", metavar='File',
                    help="File with the phenotype values for every line. One column per trait. If two phenotypes are provided (separated by commas), then the first one will be used for the Training population, and the second will be used as the Validation population. First row is header line")
parser$add_argument("-o", metavar='Directory', 
                    help = "Output directory")
parser$add_argument("-s", metavar='file', default='../geno/VEF_noMeso_annotated_repMasked_q40_s_fi_maf05_oh06_i210_imputed_rrBLUP_samples.txt',
                    help = "File with a list of samples. One sample per line. [Default %(default)s]")
parser$add_argument("-g", metavar='file', default='../geno/VEF_noMeso_annotated_repMasked_q40_s_fi_maf05_oh06_i210_imputed_rrBLUP.in',
                    help = "File with numeric genotypes. You can check \'rrBLUP\' package input for details. [Default %(default)s]")
parser$add_argument("-G", metavar='file', default='../geno/VEF_noMeso_annotated_repMasked_q40_s_fi_maf05_oh06_i210_imputed_kinship.txt',
                    help = "Kinship matrix. First row and first column are the sample names. [Default %(default)s]")

args = parser$parse_args()

setwd("/bioinfo1/projects/bean/VEF/genomic_selection/scripts")

if (any(sapply(args, is.null))){
  system('./3_run_trait_TEST.R -h')
  stop('One or more arguments are not valid. Check usage for more details.')
}

library(BGLR)

model  = unlist(strsplit(args$m,','))
traits = unlist(strsplit(args$t,','))
phen   = unlist(strsplit(args$p,','))[1]
phen2  = unlist(strsplit(args$p,','))[2]
outDir = args$o
samp = args$s
geno = args$g
Gmatrix = args$G

source("./1_getting_data_TEST.R")
source("./2_prepare_models_TEST.R")

names_list = TP_BP_partition(phen,samp,traits,phen2)

save(names_list, file = paste(outDir,'/names_list.RData',sep=''))
cat('model\ttrait\trandomPop\tcorr\tstartedAt\n')

for (i in 1:100){
  
  for (prior in model){
    
    for (trait in traits){
      
      set.seed(1234)
      combinat = names_list[[trait]]$combinat
      
      if (is.na(phen2)){
        
        if (trait %in% names(phen)){
          
          myNames = names_list[[trait]]$all_phen_gen
          
          myCorr = round( runBGLR(y = phen[myNames,],
                                  trait = trait,
                                  X = geno[myNames,],
                                  pop.split = combinat[,i],
                                  model = prior,
                                  myNames = myNames
                                  )$cor,
                          digits = 5)
          cat(prior,"\t",trait,"\tpop",i,"\t",myCorr,"\t",paste(Sys.time(),'\n'), sep = '')
        }
        
      } else {
        
        shared_traits = intersect(names(phen), names(phen2))
        
        if (trait %in% shared_traits){
          
          myNames = names_list[[trait]]$all_phen_gen2
          
          myCorr = round( runBGLR(y = phen,
                                  trait = trait,
                                  X = geno,
                                  pop.split = combinat[,i],
                                  yBP = phen2,
                                  model = prior,
                                  myNames = myNames
                                  )$cor,
                          digits = 5)
          cat(prior,"\t",trait,"\tpop",i,"\t",myCorr,"\t",paste(Sys.time(),'\n'), sep = '')
        }
      }
      }
    }
  }


