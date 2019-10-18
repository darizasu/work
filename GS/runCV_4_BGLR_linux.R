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
                    help="Comma separated list of priors to be tested. Available options are: BayesA,BayesB,BayesC,BayesRR,BLasso,FIXED,BLassof,RKHS,GBLUP")
parser$add_argument("-p", type="character", metavar='file(s)',
                    help="File with the phenotype values for every line. One column per trait. If two phenotype files are provided (separated by commas), then the first one will be used as Training population, and the second will be used as Validation population. First row is header line.")
parser$add_argument("-t", type='character', metavar='trait(s)',
                    help="Comma separated list of traits to be tested")
parser$add_argument("-o", metavar='directory', default='.',
                    help = "Output directory. [Default Current working directory]")
parser$add_argument("-g", metavar='file',
                    help = "File with numeric genotypes. You can check \'rrBLUP\' package input for details.")
parser$add_argument("-s", metavar='file',
                    help = "File with a list of samples. One sample per line.")
parser$add_argument("-G", metavar='file', default='NA',
                    help = "Kinship matrix. First row and first column are the sample names.")
parser$add_argument("-f", metavar='integer', default='30',
                    help = "Percentage of the total population to be used as validation population. This value must be between 0 - 100. [Default %(default)s]")
parser$add_argument("-I", metavar='integer', default='100',
                    help = "Number of random partition populations to be used for prediction ability assessment. Any positive integer is accepted. Proceed with caution here, the more populations the longer it takes to complete the whole analysis. [Default %(default)s]")
parser$add_argument("-O", metavar='file', default='NA',
                    help = "Name of the csv file to store the final results. Not specifying this argument avoids writing such file.")
parser$add_argument("-n", metavar='.RData', default='NA',
                    help = "NOT REQUIRED. Saved workspace with population partition matrices and genotype names. This file can be retrieved from a previous run of this script. Default behavior is to create a new file.")
parser$add_argument("-r", metavar='integer', default='NA',
                    help = "NOT REQUIRED. Keep these many randomly selected SNPs from the original matrix, the others will be filtered out. The final genotype matrix will be saved in the output directory. Default behavior is to use all available SNPs.")
parser$add_argument("-c", "--combine", action="store_true", default=FALSE,
                    help="Find traits with exactly the same set of genotypes and assign them a common random population partition for cross validation. [Default %(default)s]")
parser$add_argument("-i", "--include", action="store_true", default=FALSE,
                    help="Inlcude lines that were not present in the Training population in the prediction. Only valuable when two phenotype files are provided with the '-p' option. This option could alter the TP - VP partitions, since more lines are predicted in the validation population. [Default %(default)s]")

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
Gmatrix      = args$G ;    Gmatrix = ifelse(Gmatrix    == 'NA', NA, Gmatrix)
saveTable    = args$O ;  saveTable = ifelse(saveTable  == 'NA', NA, saveTable)
names_list   = args$n ; names_list = ifelse(names_list == 'NA', NA, names_list)
rand_SNPs    = args$r ;  rand_SNPs = ifelse(rand_SNPs  == 'NA', NA, rand_SNPs)
common       = args$combine
validPop     = as.integer(args$f)
rand_pars    = as.integer(args$I)
pnl          = args$i

source('cv4BGLR.R')

pred.results <- cv4BGLR(model = model, phen = phen, traits = traits, outDir = outDir, saveTable = saveTable,
                        geno = geno, samp = samp, Gmatrix = Gmatrix, validPop = validPop, rand_pars = rand_pars,
                        names_list = names_list, phen2 = phen2, common = common, pnl = pnl, rand_SNPs = rand_SNPs)