setwd("/bioinfo1/projects/bean/VEF/genomic_selection/scripts")

samp = read.delim('../geno/VEF_noMeso_annotated_repMasked_q40_s_fi_maf05_oh06_i210_imputed_rrBLUP_samples.txt', header = F)[,1]
geno = read.table('../geno/VEF_noMeso_annotated_repMasked_q40_s_fi_maf05_oh06_i210_imputed_rrBLUP.in', row.names = as.character(samp), header = F)
phen = read.delim(phen, row.names = 1)


# Get the number of lines with genotype and phenotype data per trait

# phen_num <-numeric(0)
# for(i in 1:ncol(phen)){
#   lines_gp = intersect( samp, rownames( phen[which(!is.na(phen[,i])),] ))
#   phen_num[i] <- length( which( !is.na( phen[lines_gp, i] )))
#   names(phen_num)[i] <- names(phen)[i]
# }

# Get the list of lines that have all variables phenotyped

all_phen_list = list()

for (name in names(phen)){
  all_phen_list[[name]] = rownames(phen)[! is.na(phen[,name])]
}

all_phen_list = Reduce(intersect, all_phen_list) # Make a large intersection with Reduce

# Get the list of samples with all phenotyped variables and genotypic data.

all_phen_gen = intersect(as.character(samp), all_phen_list)

# If there is a second dataset, get the list of lines that will be analyzed from both datasets

if (!is.na(phen2)){
  
  phen2 = read.delim(phen2, row.names = 1)
  all_phen_list2 = list()
  
  for (name in names(phen2)){
    all_phen_list2[[name]] = rownames(phen2)[! is.na(phen2[,name])]
  }
  
  # Get the list of lines that have all traits phenotyped
  all_phen_list2 = Reduce(intersect, all_phen_list2)
  
  # Get the list of samples with all phenotyped variables and genotypic data.
  
  all_phen_gen2 = intersect(as.character(samp), all_phen_list2)
  
  # Get the list of samples with pheno & geno data from both datasets
  
  all_phen_gen = intersect(all_phen_gen, all_phen_gen2)
  
}

# Generate partitions: 70% TP - 30% BP

if (is.na(combinat)){
  combinat <- matrix(0, nrow=length(all_phen_gen), ncol=100, dimnames=list(row=all_phen_gen))
  
  for( i in 1:100 ){
    # 0 = train ; 1 = test
    combinat[ sample( 1:length(all_phen_gen), length(all_phen_gen)*0.3 ) , i] = 1
  }
  
  write.table(combinat, paste(outDir,'/combinat.mtx',sep=''), sep = '\t', col.names = F, quote = F)
  
} else if (is.character(combinat)){
  
  combinat = as.matrix(read.delim(combinat,row.names = 1, header = F))
  combinat = combinat[all_phen_gen,]
  
}

# All of the tables have the same lines, which comes from all_phen_gen

if(is.data.frame(phen2)){
  all_phen_gen = all_phen_gen2 # Those lines that were not phenotyped in the 1st dataset will be predicted as well
  phen2 = phen2[all_phen_gen,]
}

geno = geno[all_phen_gen,]
phen = phen[all_phen_gen,]

