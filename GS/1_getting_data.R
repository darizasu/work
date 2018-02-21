# setwd("/bioinfo1/projects/bean/VEF/genomic_selection/scripts")

samp = read.delim(samp, header = F)[,1]
geno = read.table(geno, row.names = as.character(samp), header = F)
phen = read.delim(phen, row.names = 1)

phen2 = if (!is.na(phen2)) read.delim(phen2, row.names = 1) else NA

if (rand_SNPs != 'NA'){
  
  if (ncol(geno) < as.integer(rand_SNPs)) stop('The number of SNPs to select is higher than the actual number of SNPs available')
  
  geno = geno[,sample(1:ncol(geno), as.integer(rand_SNPs))]
  write.table(geno, paste(outDir,'/Geno_matrix_',rand_SNPs,'SNPs_rrBLUP.in',sep=''),sep = ' ', col.names=F, quote=F, row.names=F)
  cat(rand_SNPs," SNPs have been selected from the original genotype matrix. The selected SNPs were saved at ", paste(outDir,'/Geno_matrix_',rand_SNPs,'SNPs_rrBLUP.in',sep=''),'\n\n')
  
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

      namesList[[trait]]$all_phen_gen = all_phen_gen
      
      # Get the size for training and validation populations
      
      namesList[[trait]]$len = c('apg1' = length(all_phen_gen),
                                 'apg2' = length(all_phen_gen2) - length(all_phen_gen))
      
      cat('Trait\tTP-len\tVP-len\n')
      cat(trait,'\t',length(all_phen_gen),'\t',length(all_phen_gen2) - length(all_phen_gen),'\n\n')
      
    } else {
      
      # Get the total size of population
      
      cat('Trait\tLines\n')
      cat(trait,'\t',length(all_phen_gen),'\n\n')
      namesList[[trait]]$len = c('apg1' = length(all_phen_gen))
      
    }
    
    namesList[[trait]]$combinat = matrix(0, nrow=length(all_phen_gen), ncol=rand_pars, dimnames=list(row=all_phen_gen))
    
    for( i in 1:rand_pars ){
      # 0 = train ; 1 = test
      namesList[[trait]]$combinat[ sample( 1:length(all_phen_gen), length(all_phen_gen) * validPop ) , i] = 1
    }
     
  }
  
  return(namesList)

}
