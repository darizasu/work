#!/home/dariza/bin/Rscript
rm(list=ls())


LoP = c("argparse","BGLR")

new.pckgs = suppressMessages( sapply(LoP, require, character.only = TRUE) )

if (sum(!new.pckgs)){
  warning("This script will try to install the packages - ",
          toString(new.pckgs), 
          " - and all its dependencies. Otherwise abort this execution",
          immediate. = TRUE)
  Sys.sleep(10)
  install.packages(pkgs=new.pckgs,
                   repos="http://cran.r-project.org", dependencies = T)
}

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
parser$add_argument("-r", metavar='integer', default='NA',
          help = "NOT REQUIRED. Keep these many randomly selected SNPs from the original matrix, the others will be filtered out. The final genotype matrix will be saved in the output directory. Default behavior is to use all available SNPs.")
parser$add_argument("-c", "--combine", action="store_true", default=FALSE,
          help="Find traits with exactly the same set of genotypes and assign them a common random population partition for cross validation. [Default %(default)s]")
parser$add_argument("-u", "--unique", action="store_true", default=FALSE,
          help="In case a second phenotype file is used as validation population, make a unique prediction between lines in the training population and the lines that are unique to the validation population. [Default %(default)s]")

args = parser$parse_args()

if (any(sapply(args, is.null))){
  
  # Display the help message and exit
  system( paste( paste(commandArgs()[1:4], collapse = " "), 
                 '--args --help', collapse = "" ) )
  stop('One or more arguments are not valid. Check usage for more details.')
}

ca = commandArgs()
ca = sapply( split( ca[6:length(ca)],
                    ceiling( seq_along(ca[6:length(ca)]) / 2) ),
             paste, collapse = " " )

message('\nCommand line arguments:\n',
        paste0( c(commandArgs()[1:5], ca), collapse = "\n   "),'\n\n')

model        = unlist(strsplit(args$m,','))
traits       = unlist(strsplit(args$t,','))
phen         = unlist(strsplit(args$p,','))[1]
phen2        = unlist(strsplit(args$p,','))[2]
outDir       = args$o
samp         = args$s
geno         = args$g
Gmatrix      = args$G
names_list   = args$n
common       = args$combine
rand_SNPs    = args$r
validPop     = as.integer(args$f)
rand_pars    = as.integer(args$i)
unq          = args$u

rand_pars = ifelse(unq, 1, rand_pars)

if (any(validPop >= 100 | validPop < 0)){
  
  stop('The validation population partition should be a value between 0 - 100.')
  
} else {
  
  validPop = validPop / 100
  
}

samp = read.delim(samp, header = F, stringsAsFactors = F)[,1]
geno = read.table(geno, row.names = as.character(samp), header = F)
phen = read.delim(phen, row.names = 1, check.names = F)

phen2 = if (!is.na(phen2)) read.delim(phen2, row.names = 1, check.names = F) else NA

if (Gmatrix == 'NA'){

  warning("Using the function 'A.mat' from the package 'rrBLUP' to calculate the kinship matrix.\n")
  source('https://raw.githubusercontent.com/cran/rrBLUP/master/R/A.mat.R')
  G = A.mat(geno)

} else {

  message("Using the file ", Gmatrix," as the kinship matrix for the population.\n")
  G = read.delim(Gmatrix, row.names = 1, header = F)
  colnames(G) = rownames(G)
}

if (rand_SNPs != 'NA'){
  
  if (ncol(geno) < as.integer(rand_SNPs)){
  stop('The number of SNPs to select is higher than the actual number of SNPs available')
  }
  
  geno = geno[,sample(1:ncol(geno), as.integer(rand_SNPs))]
  
  write.table(geno, paste0(outDir,'/Geno_matrix_',rand_SNPs,'SNPs_rrBLUP.in'),
        sep = ' ', col.names = F, quote = F, row.names = F)
  
  message(rand_SNPs," SNPs have been selected from the original genotype matrix.",
      " The selected SNPs were saved at ",
      paste0(outDir,'/Geno_matrix_',rand_SNPs,'SNPs_rrBLUP.in'),'\n')
}

source("https://raw.githubusercontent.com/darizasu/work/master/GS/TP_BP_partition.R")
source("https://raw.githubusercontent.com/darizasu/work/master/GS/runBGLR.R")

if (names_list == 'NA'){

  message('Creating the population partition matrices and subsetting genotypes with phenotypic and genotypic data available.\n')

  names_list = TP_BP_partition(myPhen = phen, myGen = samp, traits = traits, 
                 myPhen2 = phen2, validation = validPop, unq = unq)
  
  save(names_list, file = paste0(outDir,'/names_list.RData'))
  message('The population partition matrices and genotype names have been saved at:',
          paste0(outDir,'/names_list.RData'),'\n')
  
} else {
  
  namesListFile = names_list # The string containing the location of the .Rdata file is saved in namesListFile
  message('Loading population partition matrices and genotype names from:',names_list,'\n')
  load(names_list) # The 'name_list' object is no longer a string with the location of the .Rdata file, instead is the list that was contained in the file.
  
    for (trait in names(names_list)){
  
    if (ncol(names_list[[trait]]$combinat) < rand_pars){
      stop(paste('There are only',ncol(names_list[[trait]]$combinat),
                 'population partitions available at',namesListFile,
                 'for the trait',trait,'and you requested',rand_pars,
                 "partitions with '-i' option."))
    }
  }
}

if (common){
  
  rand_pops = lapply(names_list, `[[`, 'l2bu')
  unique_pops = unique(rand_pops)
  
  pops_n = sapply(rand_pops, function(x) which( sapply(unique_pops, identical, x) ))
  
  message('\n\nA common population partition matrix will be assigned to traits with exactly the same set of genotypes:\n')

  for (i in unique(pops_n)){

    pop = names(pops_n[pops_n == i])
    first = names_list[[ pop[1] ]]$combinat
    message(paste0(' ├─ ', pop, coll = '\n'))

    for (j in pop) names_list[[ j ]]$combinat = first
  }
}

sbt = 1
message('\n\nGenotype composition for each trait:\n')

for (trait in traits){

  if (is.data.frame(phen2)){
  
    if (sbt) message('  ├─\tTrait\t=\tTP-length\tVP-length\n')
    message(' ── ',trait,'\t=\t', length(names_list[[trait]]$linePG),'\t=\t',
            length(names_list[[trait]]$linePG_vp))
    sbt = 0
  
  } else {
  
    if (sbt) message('  ├─\tTrait\t=\tLines\n')
    message('  ──  ',trait,'\t=\t',length(names_list[[trait]]$linePG))
    sbt = 0
  }
}

cat('model\ttrait\trandomPop\tcorr\tfinishedAt\n')

myTable = data.frame(model = NA, trait = NA, pred = NA, obs = NA)

for (i in 1:rand_pars){
  
  for (prior in model){
  
    for (trait in traits){
    
      set.seed(1234)
      combinat = names_list[[trait]]$combinat
      myNames  = names_list[[trait]]$l2bu

      if (! is.data.frame(phen2)){
    
        if (trait %in% names(phen)){
      
          myCorr = round( runBGLR(y = phen[myNames,],
                                  trait = trait,
                                  X = geno[myNames,],
                                  pop.split = combinat[,i],
                                  model = prior,
                                  myNames = myNames,
                                  G = G[myNames,myNames]
                                  )$cor,
                          digits = 5)
          cat(prior, trait, paste0('pop',i), myCorr, paste0(Sys.time(),'\n'))
        }
    
      } else {
    
        shared_traits = intersect(names(phen), names(phen2))
    
        if (trait %in% shared_traits){
      
          myCorr = runBGLR(y = phen,
                           trait = trait,
                           X = geno,
                           pop.split = combinat[,i],
                           yObs = phen2,
                           model = prior,
                           myNames = myNames,
                           G = G[myNames,myNames]
                           )$cor
          cat(prior, trait, paste0('pop',i), myCorr, paste0(Sys.time(),'\n'))
          # myTable = rbind(myTable, myCorr)
        }
      }
    }
  }
}

# write.csv(myTable, 'myTable.csv', quote = F, row.names = F)
