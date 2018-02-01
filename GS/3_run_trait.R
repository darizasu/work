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

parser$add_argument("-m", type='character', metavar='model(s)',
                    help="Comma separated list of priors to be tested. Available options are BayesA,BayesB,BayesC,BayesRR,BLasso,FIXED,BLassof,RKHS,GBLUP")
parser$add_argument("-t", type='character', metavar='trait(s)',
                    help="Comma separated list of traits to be tested")
parser$add_argument("-p", type="character", metavar='file',
                    help="File with the phenotype values for every line. One column per trait. If two phenotypes are provided (separated by commas), then the first one will be used for the Training population, and the second will be used as the Validation population. First row is header line")
parser$add_argument("-o", metavar='directory', default='.',
                    help = "Output directory. [Default /bioinfo1/projects/bean/VEF/genomic_selection/scripts]")
parser$add_argument("-s", metavar='file', default='../geno/VEF_noMeso_annotated_repMasked_q40_s_fi_oh06_i210_maf0.05_imputed_rrBLUP_samples.txt',
                    help = "File with a list of samples. One sample per line. [Default %(default)s]")
parser$add_argument("-g", metavar='file', default='../geno/VEF_noMeso_annotated_repMasked_q40_s_fi_oh06_i210_maf0.05_imputed_rrBLUP.in',
                    help = "File with numeric genotypes. You can check \'rrBLUP\' package input for details. [Default %(default)s]")
parser$add_argument("-G", metavar='file', default='../geno/VEF_noMeso_annotated_repMasked_q40_s_fi_oh06_i210_maf0.05_imputed_kinship.txt',
                    help = "Kinship matrix. First row and first column are the sample names. [Default %(default)s]")
parser$add_argument("-n", metavar='.RData', default='NA',
                    help = "NOT REQUIRED. Saved workspace with population partition matrices and genotype names. This file can be retrieved from a previous run of this script. Default behavior is to create a new file.")
parser$add_argument("-r", metavar='character', default='NA',
                    help = "NOT REQUIRED. Keep these many randomly selected SNPs from the original matrix, the others will be filtered out. The final genotype matrix will be saved in the output directory. Default behavior is to use all available SNPs.")
parser$add_argument("-f", metavar='integer', default='70',
                    help = "NOT REQUIRED. Percentage of the total population to be used as training population. This value should be between 0 - 100. [Default %(default)s]")

args = parser$parse_args()


if (any(sapply(args, is.null))){
  
  system('/bioinfo1/projects/bean/VEF/genomic_selection/scripts/3_run_trait.R -h')
  stop('One or more arguments are not valid. Check usage for more details.')
  
}

cat('\nCommand line arguments:\n',commandArgs(),'\n\n')

setwd("/bioinfo1/projects/bean/VEF/genomic_selection/scripts")

library(BGLR)

model  = unlist(strsplit(args$m,','))
traits = unlist(strsplit(args$t,','))
phen   = unlist(strsplit(args$p,','))[1]
phen2  = unlist(strsplit(args$p,','))[2]
outDir = args$o
samp = args$s
geno = args$g
Gmatrix = args$G
names_list = args$n
rand_SNPs = args$r
trainPop = as.integer(args$f)

if (any(trainPop > 100 | trainPop <= 0)){
  stop('The training population partition should be a value between 0 - 100.')
} else {
  trainPop = trainPop / 100
}

source("/bioinfo1/projects/bean/VEF/genomic_selection/scripts/1_getting_data.R")
source("/bioinfo1/projects/bean/VEF/genomic_selection/scripts/2_prepare_models.R")

if (names_list == 'NA'){
  
  names_list = TP_BP_partition(phen,samp,traits,phen2)
  save(names_list, file = paste(outDir,'/names_list.RData',sep=''))
  cat('The population partition matrices and genotype names have been saved in:',
      paste(outDir,'/names_list.RData',sep=''),'\n\n')
  
} else {
  
  cat('Loading population partition matrices and genotype names from:', names_list,'\n\n')
  load(names_list)
  
}

cat('model\ttrait\trandomPop\tcorr\tfinishedAt\n')

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
                                  myNames = myNames,
                                  G = Gmatrix
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
                                  myNames = myNames,
                                  G = Gmatrix
                                  )$cor,
                          digits = 5)
          cat(prior,"\t",trait,"\tpop",i,"\t",myCorr,"\t",paste(Sys.time(),'\n'), sep = '')
        }
      }
      }
    }
  }


