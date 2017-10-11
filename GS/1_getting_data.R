setwd("/bioinfo1/projects/bean/VEF/genomic_selection/scripts")

samp = read.delim('../geno/VEF_noMeso_annotated_repMasked_q40_s_fi_maf05_oh06_i210_imputed_rrBLUP_samples.txt', header = F)[,1]
geno = read.table('../geno/VEF_noMeso_annotated_repMasked_q40_s_fi_maf05_oh06_i210_imputed_rrBLUP.in', row.names = as.character(samp), header = F)
phen = read.delim(phen, row.names = 1)
if (!is.na(phen2)){
  phen2 = read.delim(phen2, row.names = 1)
  combinat = NA
}


# Get the number of lines with genotype and phenotype data per trait

# phen_num <-numeric(0)
# for(i in 1:ncol(phen)){
#   lines_gp = intersect( samp, rownames( phen[which(!is.na(phen[,i])),] ))
#   phen_num[i] <- length( which( !is.na( phen[lines_gp, i] )))
#   names(phen_num)[i] <- names(phen)[i]
# }

# Get the list of lines that have all phenotyped variables

all_phen_list = list()

for (name in names(phen)){
  all_phen_list[[name]] = rownames(phen)[! is.na(phen[,name])]
}

all_phen_list = Reduce(intersect, all_phen_list)

# Get the list of samples with all phenotyped variables and genotypic data.

all_phen_gen = intersect(as.character(samp), all_phen_list)

# Generate partitions: 70% TP - 30% BP

if (is.na(phen2) && is.na(combinat)){
  combinat <- matrix(0, nrow=length(all_phen_gen), ncol=100, dimnames=list(row=all_phen_gen))
  
  for( i in 1:100 ){
    # 0 = train ; 1 = test
    combinat[ sample( 1:length(all_phen_gen), length(all_phen_gen)*0.3 ) , i] = 1
  }
  
  # hist(rowSums(combinat[,1:100]),100)
  write.table(combinat, paste(outDir,'/combinat.mtx',sep=''), sep = '\t', col.names = F, quote = F)
  
} else if (is.character(combinat)){
  
  combinat = as.matrix(read.delim(combinat,row.names = 1, header = F))
  
}

geno = geno[all_phen_gen,]
phen = phen[all_phen_gen,]
if (is.matrix(combinat)) combinat = combinat[all_phen_gen,]

# X = geno
# Z = scale(X, center = T, scale = T)
# G = tcrossprod(Z) / ncol(Z)
# G = read.delim('VEF_annotated_repMasked_q40_s_fi_maf05_oh06_i290_kinship.txt', row.names = 1)
# G = as.matrix(G[all_phen_gen,all_phen_gen])

trait.cors = data.frame(BayesA = NA, BayesB = NA, BayesC = NA, BayesRR = NA, BLasso = NA, BLassof = NA, RKHS = NA, GBLUP = NA, FIXED = NA)

cors = list()

for (trait in traits){
  cors[[trait]] = trait.cors
}

