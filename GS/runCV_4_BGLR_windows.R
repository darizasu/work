rm(list=ls())

# model  = Priors to be tested in a vector of characters. Available options are BayesA,BayesB,BayesC,BayesRR,BLasso,FIXED,BLassof,RKHS,GBLUP
# model  = c('BayesA','BayesB','BayesC','BayesRR','BLasso','FIXED','BLassof','RKHS','GBLUP')
model  = c('RKHS','GBLUP')

# traits = Traits to run genomic prediction, in a vector of characters. The names must be the same than the ones in the following 'phen' file.
traits = c('wholepop_DF','Fe','Zn','SdWt','subsamp_DF','wac','p80CKT','p50CKT','p92CKT','p100CKT','POM','TSW','100SW','YdPl_ha2','subsamp_POM','subsamp_TSW','subsamp_100SW','subsamp_YdPl_ha2','raw_p80CKT')

# phen   = File with phenotypic data. First column is for genotype names. Then, one column per trait
phen   = 'GS_7.txt'

# phen2  = NOT REQUIRED. File with phenotypic data from a different year, to carry out year to year prediction.
phen2  = NA

# outDir = Output directory
outDir = getwd()

# samp = File with the list of genotype names present in the genotypic matrix. One per line.
samp = 'MII_repMasked_q40_s_fi_maf05_oh06_i130_noScaffolds_imputed_rrBLUP_samples.txt'

# geno = A matrix with genotypic data coded as (-1,0,1) or (0,1,2) for homozygous reference allele, heterozygous, and homozygous alternative allele.
# One column per marker and one row per genotype. No headers. The order of genotypes is the same order than that of the previous 'samp' file
geno = 'MII_repMasked_q40_s_fi_maf05_oh06_i130_noScaffolds_imputed_rrBLUP.in'

# Gmatrix = Kinship matrix. No headers. First column is for the genotype names
Gmatrix = 'MII_repMasked_q40_s_fi_maf05_oh06_i130_noScaffolds_kinship.txt'

# validPop = Percentage of the total population to be used as validation population. This value must be between 0 - 100.
validPop = 30

# rand_pars = Number of random partition populations to be used for prediction ability assessment. Any positive integer is accepted.
# Proceed with caution here, the more populations the longer it takes to complete the whole analysis.
rand_pars = 50

# names_list = NOT REQUIRED. Saved workspace with population partition matrices and genotype names. 
# This file can be retrieved from a previous run of this script. Default behavior is to create a new file.
names_list = NA

# Find traits with exactly the same set of genotypes and assign them a common random population partition for cross validation.
common = TRUE

# rand_SNPs = NOT REQUIRED. Keep these many randomly selected SNPs from the original matrix, the others will be filtered out.
# The final genotype matrix will be saved in the output directory. Default behavior is to use all available SNPs.
rand_SNPs = NA


#### ---------------------------------------------------------------------- ####
#### ----------- DO NOT MODIFY ANYTHING FROM THIS POINT FORWARD ----------- ####
#### ---------------------------------------------------------------------- ####

list.of.packages = "BGLR"
new.pckgs = list.of.packages[!(list.of.packages %in% rownames(installed.packages()))]

if (length(new.pckgs)){
  cat("\nWarning message:\nThis script will try to install the packages",new.pckgs,"and all its dependencies. Otherwise abort this execution\n")
  Sys.sleep(10)
  install.packages(pkgs=new.pckgs, repos="http://cran.r-project.org", dependencies = T)
}

library(BGLR)

if (any(validPop > 100 | validPop <= 0)){
  
  stop('The validation population partition should be a value between 0 - 100.')
  
} else {
  
  validPop = validPop / 100
  
}

samp = read.delim(samp, header = F)[,1]
geno = read.table(geno, row.names = as.character(samp), header = F)
phen = read.delim(phen, row.names = 1)

phen2 = if (!is.na(phen2)) read.delim(phen2, row.names = 1) else NA

if (!is.na(rand_SNPs)){
  
  if (ncol(geno) < as.integer(rand_SNPs)) stop('The number of SNPs to select is higher than the actual number of SNPs available')
  
  geno = geno[,sample(1:ncol(geno), as.integer(rand_SNPs))]
  write.table(geno, paste(outDir,'/Geno_matrix_',rand_SNPs,'SNPs_rrBLUP.in',sep=''),sep = ' ', col.names=F, quote=F, row.names=F)
  cat(rand_SNPs," SNPs have been selected from the original genotype matrix. The selected SNPs were saved at ", paste(outDir,'/Geno_matrix_',rand_SNPs,'SNPs_rrBLUP.in',sep=''),'\n\n')
  
}

source("https://raw.githubusercontent.com/darizasu/work/master/GS/TP_BP_partition.R")
source("https://raw.githubusercontent.com/darizasu/work/master/GS/runBGLR.R")

if (is.na(names_list)){
  
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
    
    cat('A common population partition matrix was assigned to traits with exactly the same set of genotypes:\n')
    print(lapply(pops_n, names),end='\n\n')
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
predResults = data.frame(model=character(),
                         trait=character(),
                         randomPop=character(),
                         corr=numeric(),
                         finishedAt=character())

for (i in 1:rand_pars){
  
  pop = paste0('pop',i)
  
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
          tmp_DF = data.frame(prior,trait,pop,myCorr,Sys.time())
          names(tmp_DF) = names(predResults)
          predResults = rbind(predResults, tmp_DF)
          
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
          tmp_DF = data.frame(prior,trait,pop,myCorr,Sys.time())
          names(tmp_DF) = names(predResults)
          predResults = rbind(predResults, tmp_DF)
          
        }
      }
    }
  }
}

write.csv(predResults, 'Prediction_results.csv', quote=F, row.names=F)

