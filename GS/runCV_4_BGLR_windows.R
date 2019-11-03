

# model  <- Priors to be tested in a vector of characters. Available options are: BayesA,BayesB,BayesC,BayesRR,BLasso,FIXED,BLassof,RKHS,GBLUP
# model  <- c('BayesA','BayesB','BayesC','BayesRR','BLasso','FIXED','BLassof','RKHS','GBLUP')
model  <- c('RKHS')

# phen   <- File with phenotypic data. First column is for genotype names. Then, one column per trait
phen   <- 'D:/OneDrive - CGIAR/Daniel/Winny_Amongi/GS Data MEANS.txt'

# traits <- Traits to run genomic prediction, in a vector of characters. The names must be the same than the ones in the 'phen' file.
# traits <- c('ALSF','ALSFP','ANTFL','ANTFP','CBBFL','CBBFP','RUSTFL','RUSTFP','PLNTVIG','DF','DPM','YDHA','FESEED','ZNSEED')
traits <- c('YDHA','FESEED','ZNSEED')

# out_dir <- Output directory
out_dir <- 'D:/OneDrive - CGIAR/Daniel/Winny_Amongi'

# Name of the csv file to store the final results. NA avoids writing such file.
save_table <- 'D:/OneDrive - CGIAR/Daniel/Winny_Amongi/Prediction_results.csv'

# geno <- A matrix with genotypic data coded as (-1,0,1) or (0,1,2) for homozygous reference allele, heterozygous, and homozygous alternative allele.
# One column per marker and one row per genotype. No headers. The order of genotypes is the same order than that of the next 'samp' file
geno <- 'D:/OneDrive - CGIAR/Daniel/Winny_Amongi/GBS_Africa_population_repMasked_mi300_annotated_maf0.1_BeagleImpute_annotated_rrBLUP.in'

# samp <- File with the list of genotype names present in the genotypic matrix. One per line.
samp <- 'D:/OneDrive - CGIAR/Daniel/Winny_Amongi/GBS_Africa_population_repMasked_mi300_annotated_maf0.1_BeagleImpute_annotated_rrBLUP_samples.txt'

# Gmatrix <- Kinship matrix. No headers. First column is for the genotype names. In case no Gmatrix is provided, the function 'A.mat' from the package rrBLUP is used to calculate it.
Gmatrix <- NA

# validation_perc <- Percentage of the total population to be used as validation population. This value must be an integer between 0 - 100.
validation_perc <- 30

# rand_pars <- Number of random partition populations to be used for prediction ability assessment. Any positive integer is accepted.
# Proceed with caution here, the more populations the longer it takes to complete the whole analysis.
rand_pars <- 2

# names_list <- NOT REQUIRED. Saved workspace with population partition matrices and genotype names. 
# This file can be retrieved from a previous run of this script. Default behavior is to create a new file.
names_list <- NA

# phen2  <- NOT REQUIRED. File with phenotypic data from a different trial. This file will be used as the validation set, and 'phen' will be used as the training set.
phen2  <- NA

# Find traits with exactly the same set of genotypes and assign them a common random population partition for cross validation.
common <- TRUE

# Inlcude lines that were not present in the Training population in the prediction.
# Only valuable when two phenotype files are provided with the '-p' option.
# This option could alter the TP - VP partitions, since more lines are predicted in the validation population.
pnl <- FALSE

# rand_SNPs <- NOT REQUIRED. Keep these many randomly selected SNPs from the original matrix, the others will be filtered out.
# The final genotype matrix will be saved in the output directory. If NA, default behavior is to use all available SNPs.
rand_SNPs <- NA

source('https://raw.githubusercontent.com/darizasu/work/master/GS/cv4BGLR.R')

pred_results <- cv4BGLR(model = model, phen = phen, traits = traits, out_dir = out_dir, save_table = save_table,
                        geno = geno, samp = samp, Gmatrix = Gmatrix, validation_perc = validation_perc, rand_pars = rand_pars,
                        names_list = names_list, phen2 = phen2, common = common, pnl = pnl, rand_SNPs = rand_SNPs)
