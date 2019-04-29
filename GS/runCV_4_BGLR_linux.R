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
                    help = "Output directory. [Default Current working directory]")
parser$add_argument("-s", metavar='file',
                    help = "File with a list of samples. One sample per line.")
parser$add_argument("-g", metavar='file',
                    help = "File with numeric genotypes. You can check \'rrBLUP\' package input for details.")
parser$add_argument("-G", metavar='file', default='NA',
                    help = "Kinship matrix. First row and first column are the sample names.")
parser$add_argument("-f", metavar='integer', default='30',
                    help = "Percentage of the total population to be used as validation population. This value must be between 0 - 100. [Default %(default)s]")
parser$add_argument("-i", metavar='integer', default='100',
                    help = "Number of random partition populations to be used for prediction ability assessment. Any positive integer is accepted. Proceed with caution here, the more populations the longer it takes to complete the whole analysis. [Default %(default)s]")
parser$add_argument("-n", metavar='.RData', default='NA',
                    help = "NOT REQUIRED. Saved workspace with population partition matrices and genotype names. This file can be retrieved from a previous run of this script. Default behavior is to create a new file.")
parser$add_argument("-c", "--combine", action="store_true", default=FALSE,
                    help="Find traits with exactly the same set of genotypes and assign them a common random population partition for cross validation. [Default %(default)s]")
parser$add_argument("-r", metavar='integer', default='NA',
                    help = "NOT REQUIRED. Keep these many randomly selected SNPs from the original matrix, the others will be filtered out. The final genotype matrix will be saved in the output directory. Default behavior is to use all available SNPs.")

args = parser$parse_args()

if (any(sapply(args, is.null))){
  
  # Display the help message and exit
  system( paste( paste(commandArgs()[1:4], collapse = " "), '--args --help', collapse = "" ) )
  stop('One or more arguments are not valid. Check usage for more details.')
  
}

ca = commandArgs()
ca = sapply( split( ca[6:length(ca)],
                    ceiling( seq_along(ca[6:length(ca)]) / 2) ),
             paste, collapse = " " )

cat('\nCommand line arguments:\n',
    paste0( c(commandArgs()[1:5], ca), collapse = "\n   "),'\n\n')

library(BGLR)

model      = unlist(strsplit(args$m,','))
traits     = unlist(strsplit(args$t,','))
phen       = unlist(strsplit(args$p,','))[1]
phen2      = unlist(strsplit(args$p,','))[2]
outDir     = args$o
samp       = args$s
geno       = args$g
Gmatrix    = args$G
names_list = args$n
common     = args$combine
rand_SNPs  = args$r
validPop   = as.integer(args$f)
rand_pars  = as.integer(args$i)

if (any(validPop > 100 | validPop <= 0)){
  
  stop('The validation population partition should be a value between 0 - 100.')
  
} else {
  
  validPop = validPop / 100
  
}

samp = read.delim(samp, header = F, stringsAsFactors = F)[,1]
geno = read.table(geno, row.names = as.character(samp), header = F)
phen = read.delim(phen, row.names = 1, check.names = F)

phen2 = if (!is.na(phen2)) read.delim(phen2, row.names = 1, check.names = F) else NA

if (Gmatrix == 'NA'){
  
  cat("\nUsing the function 'A.mat' from the package 'rrBLUP' to calculate the kinship matrix.\n\n")
  source('https://raw.githubusercontent.com/cran/rrBLUP/master/R/A.mat.R')
  G = A.mat(geno)
  
} else {
  
  cat("\nUsing the file ", Gmatrix," as the kinship matrix for the population.\n\n")
  G = read.delim(Gmatrix, row.names = 1, header = F)
  colnames(G) = rownames(G)
}

if (rand_SNPs != 'NA'){
  
  if (ncol(geno) < as.integer(rand_SNPs)) stop('The number of SNPs to select is higher than the actual number of SNPs available')
  
  geno = geno[,sample(1:ncol(geno), as.integer(rand_SNPs))]
  
  write.table(geno, paste(outDir,'/Geno_matrix_',rand_SNPs,'SNPs_rrBLUP.in',sep=''),sep = ' ',
              col.names=F, quote=F, row.names=F)
  
  cat(rand_SNPs," SNPs have been selected from the original genotype matrix. The selected SNPs were saved at ",
      paste(outDir,'/Geno_matrix_',rand_SNPs,'SNPs_rrBLUP.in',sep=''),'\n\n')
  
}

source("https://raw.githubusercontent.com/darizasu/work/master/GS/TP_BP_partition.R")
source("https://raw.githubusercontent.com/darizasu/work/master/GS/runBGLR.R")

if (names_list == 'NA'){
  
  names_list = TP_BP_partition(phen,samp,traits,phen2)
  
  if (common){
    
    rand_pops = lapply(names_list, `[[`, 'all_phen_gen')
    unique_pops = unique(rand_pops)
    
    pops_n = sapply(unique_pops, function(x) which( sapply(rand_pops, identical, x) ))
    
    for (i in 1:length(pops_n)){
      
      first = names_list[[ pops_n[[i]][1] ]]$combinat
      
      for (j in pops_n[[i]]){
        
        names_list[[j]]$combinat = first
        
      }
    }
    
    cat('\nA common population partition matrix was assigned to traits with exactly the same set of genotypes:\n')
    print(lapply(pops_n, names), end='\n\n')
  }
  
  save(names_list, file = paste(outDir,'/names_list.RData',sep=''))
  cat('The population partition matrices and genotype names have been saved at:',
      paste(outDir,'/names_list.RData',sep=''),'\n\n')
  
} else {
  
  namesListFile = names_list # The string containing the location of the .Rdata file is saved in namesListFile
  cat('Loading population partition matrices and genotype names from:',names_list,'\n\n')
  load(names_list) # The 'name_list' object is no longer a string with the location of the .Rdata file, instead is the list that was contained in the file.
  
  for (trait in names(names_list)){
    
    if (ncol(names_list[[trait]]$combinat) < rand_pars){
      stop(paste('There are only',ncol(names_list[[trait]]$combinat),'population partitions available at',namesListFile,'for the trait',trait,'and you requested',rand_pars, "partitions with '-i' option."))
    }
    
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

for (i in 1:rand_pars){
  
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
                                  G = G[myNames,myNames]
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
                                  G = G[myNames,myNames]
          )$cor,
          digits = 5)
          cat(prior,"\t",trait,"\tpop",i,"\t",myCorr,"\t",paste(Sys.time(),'\n'), sep = '')
        }
      }
    }
  }
}


