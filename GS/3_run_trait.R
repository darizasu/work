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
parser$add_argument("-p", type="character", metavar='file(s)',
                    help="File with the phenotype values for every line. One column per trait. If two phenotype files are provided (separated by commas), then the first one will be used as Training population, and the second will be used as Validation population. First row is header line.")
parser$add_argument("-o", metavar='directory', default='.',
                    help = "Output directory. [Default /bioinfo1/projects/bean/VEF/genomic_selection/scripts]")
parser$add_argument("-s", metavar='file', default='../geno/VEF_noMeso_annotated_repMasked_q40_s_fi_oh06_i210_maf0.05_imputed_rrBLUP_samples.txt',
                    help = "File with a list of samples. One sample per line. [Default %(default)s]")
parser$add_argument("-g", metavar='file', default='../geno/VEF_noMeso_annotated_repMasked_q40_s_fi_oh06_i210_maf0.05_imputed_rrBLUP.in',
                    help = "File with numeric genotypes. You can check \'rrBLUP\' package input for details. [Default %(default)s]")
parser$add_argument("-G", metavar='file', default='../geno/VEF_noMeso_annotated_repMasked_q40_s_fi_oh06_i210_maf0.05_imputed_kinship.txt',
                    help = "Kinship matrix. First row and first column are the sample names. [Default %(default)s]")
parser$add_argument("-f", metavar='integer', default='30',
                    help = "Percentage of the total population to be used as validation population. This value must be between 0 - 100. [Default %(default)s]")
parser$add_argument("-n", metavar='.RData', default='NA',
                    help = "NOT REQUIRED. Saved workspace with population partition matrices and genotype names. This file can be retrieved from a previous run of this script. Default behavior is to create a new file.")
parser$add_argument("-r", metavar='integer', default='NA',
                    help = "NOT REQUIRED. Keep these many randomly selected SNPs from the original matrix, the others will be filtered out. The final genotype matrix will be saved in the output directory. Default behavior is to use all available SNPs.")


args = parser$parse_args()


if (any(sapply(args, is.null))){
  
  system('/bioinfo1/projects/bean/VEF/genomic_selection/scripts/3_run_trait.R -h') # This script should be run no matter where it is located. Change this.
  stop('One or more arguments are not valid. Check usage for more details.')
  
}

cat('\nCommand line arguments:\n',commandArgs(),'\n\n')

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
validPop = as.integer(args$f)

if (any(validPop > 100 | validPop <= 0)){
  
  stop('The validation population partition should be a value between 0 - 100.')
  
} else {
  
  validPop = validPop / 100
  
}

source("/bioinfo1/projects/bean/VEF/genomic_selection/scripts/1_getting_data.R")   # This script should be run no matter where it is located. Change this.
source("/bioinfo1/projects/bean/VEF/genomic_selection/scripts/2_prepare_models.R") # This script should be run no matter where it is located. Change this.

if (names_list == 'NA'){
  
  names_list = TP_BP_partition(phen,samp,traits,phen2)
  save(names_list, file = paste(outDir,'/names_list.RData',sep=''))
  cat('The population partition matrices and genotype names have been saved at:',
      paste(outDir,'/names_list.RData',sep=''),'\n\n')
  
} else {
  
  cat('Loading population partition matrices and genotype names from:',names_list,'\n\n')
  load(names_list)
  
  for (trait in names(names_list)){
    
    if (is.data.frame(phen2)){
      
      cat('Trait\tTP-len\tVP-len\n')
      cat(trait,'\t',names_list[[trait]]$len['apg1'],'\t',names_list[[trait]]$len['apg2'],'\n\n')
      
    } else {
      
      cat('Trait\tLines\n')
      cat(trait,'\t',names_list[[trait]]$len['apg1'],'\n\n')
      
    }
    
  }
  
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


