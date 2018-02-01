setwd("/bioinfo1/projects/bean/VEF/genomic_selection/scripts")

samp = read.delim(samp, header = F)[,1]
geno = read.table(geno, row.names = as.character(samp), header = F)
phen = read.delim(phen, row.names = 1)

phen2 = if (!is.na(phen2)) read.delim(phen2, row.names = 1) else NA

if (rand_SNPs != 'NA'){
  
  if (ncol(geno) < as.integer(rand_SNPs)) stop('The number of SNPs to select is higher than the actual number of SNPs available')
  
  geno = geno[,sample(1:ncol(geno), as.integer(rand_SNPs))]
  write.table(geno, paste(outDir,'/Geno_matrix_',rand_SNPs,'SNPs_rrBLUP.in',sep=''),sep = ' ', col.names=F, quote=F, row.names=F)
  cat(rand_SNPs," SNPs have been selected from the original genotype matrix. The selected SNPs were saved in ", paste(outDir,'/Geno_matrix_',rand_SNPs,'SNPs_rrBLUP.in',sep=''),'\n\n')
  
}
# Get the number of lines with genotype and phenotype data per trait

# phen_num <-numeric(0)
# for(i in 1:ncol(phen)){
#   lines_gp = intersect( samp, rownames( phen[which(!is.na(phen[,i])),] ))
#   phen_num[i] <- length( which( !is.na( phen[lines_gp, i] )))
#   names(phen_num)[i] <- names(phen)[i]
# }

# A function that returns you the all_phen_gen for a given trait and its corresponding combinat

TP_BP_partition = function (myPhen,myGen,traits,myPhen2){
  
  namesList = list()
  
  # Avoid those traits that are not present in myPhen
  
  traits = traits[traits %in% names(myPhen)]
  if (is.data.frame(myPhen2)) traits = traits[traits %in% names(myPhen2)]
  
  for (trait in traits){
   
    # Get the list of lines that have all variables phenotyped
    all_phen_list = rownames(myPhen)[! is.na(myPhen[,trait])]
    
    # Get the list of samples with all phenotyped variables and genotypic data
    all_phen_gen = intersect(as.character(myGen), all_phen_list)
    
    namesList[[trait]]$all_phen_gen = all_phen_gen
    
    # If there is a 2nd dataset, get the list of lines that will be analyzed from both datasets
    if (is.data.frame(myPhen2)){
      
      
      all_phen_list2 = rownames(myPhen2)[! is.na(myPhen2[,trait])]
      
      all_phen_gen2 = intersect(as.character(samp), all_phen_list2)
      
      namesList[[trait]]$all_phen_gen2 = all_phen_gen2
      
      # Get the list of samples with pheno & geno data from both datasets
      
      all_phen_gen = intersect(all_phen_gen, all_phen_gen2)
      # cat('trait\tLines_in_both_datasets\tUnique_in_2nd_dataset\n')
      cat(trait,'\t',length(all_phen_gen),'\t',length(all_phen_gen2) - length(all_phen_gen),'\n')
      
      namesList[[trait]]$all_phen_gen = all_phen_gen
      
      
    } else {
      
      # cat('trait\tLines_with_phen&gen\n')
      cat(trait,'\t',length(all_phen_gen),'\n')
      
    }
    
    namesList[[trait]]$combinat = matrix(0, nrow=length(all_phen_gen), ncol=100, dimnames=list(row=all_phen_gen))
    
    for( i in 1:100 ){
      # 0 = train ; 1 = test
      namesList[[trait]]$combinat[ sample( 1:length(all_phen_gen), length(all_phen_gen) * validPop ) , i] = 1
    }
     
  }
  
  return(namesList)

}

# # Get the list of lines that have all variables phenotyped
# 
# all_phen_list = list()
# 
# for (name in names(phen)){
#   all_phen_list[[name]] = rownames(phen)[! is.na(phen[,name])]
# }
# 
# all_phen_list = Reduce(intersect, all_phen_list) # Make a large intersection with Reduce
# 
# # Get the list of samples with all phenotyped variables and genotypic data.
# 
# all_phen_gen = intersect(as.character(samp), all_phen_list)
# 
# # If there is a second dataset, get the list of lines that will be analyzed from both datasets
# 
# if (!is.na(phen2)){
#   
#   phen2 = read.delim(phen2, row.names = 1)
#   all_phen_list2 = list()
#   
#   for (name in names(phen2)){
#     all_phen_list2[[name]] = rownames(phen2)[! is.na(phen2[,name])]
#   }
#   
#   # Get the list of lines that have all traits phenotyped
#   all_phen_list2 = Reduce(intersect, all_phen_list2)
#   
#   # Get the list of samples with all phenotyped variables and genotypic data.
#   
#   all_phen_gen2 = intersect(as.character(samp), all_phen_list2)
#   
#   # Get the list of samples with pheno & geno data from both datasets
#   
#   all_phen_gen = intersect(all_phen_gen, all_phen_gen2)
#   
# }
# 
# # Generate partitions: 70% TP - 30% BP
# 
# if (is.na(combinat)){
#   combinat <- matrix(0, nrow=length(all_phen_gen), ncol=100, dimnames=list(row=all_phen_gen))
#   
#   for( i in 1:100 ){
#     # 0 = train ; 1 = test
#     combinat[ sample( 1:length(all_phen_gen), length(all_phen_gen)*0.3 ) , i] = 1
#   }
#   
#   write.table(combinat, paste(outDir,'/combinat.mtx',sep=''), sep = '\t', col.names = F, quote = F)
#   
# } else if (is.character(combinat)){
#   
#   combinat = as.matrix(read.delim(combinat,row.names = 1, header = F))
#   combinat = combinat[all_phen_gen,]
#   
# }
# 
# # All of the tables have the same lines, which comes from all_phen_gen
# 
# if(is.data.frame(phen2)){
#   all_phen_gen = all_phen_gen2 # Those lines that were not phenotyped in the 1st dataset will be predicted as well
#   phen2 = phen2[all_phen_gen,]
# }
# 
# geno = geno[all_phen_gen,]
# phen = phen[all_phen_gen,]

