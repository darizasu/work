

TP_BP_partition = function (myPhen,myGen,traits,myPhen2,validation,predNewLines){
  
  # Object . A function that returns the genotype IDs that have phenotypic and 
  #          genotypic information for each requested trait. Based on this information,
  #          it produces the random assignment matrices for training and validation populations.
  # Input . myPhen . 
  # Input . myGen . 
  # Input . traits . 
  # Input . myPhen2 .
  # Output . A list with as many items as valid traits there are. 
  #          Each trait is a sublist containing a vector with genotype IDs (all_phen_gen), 
  #          the length of that vector (len) and the Training and Validation partition matrices (combinat).
  # Authors: dariza and jdelahoz
  #   Last modified: November 2, 2018
  
  namesList = list()
  
  # Avoid those traits that are not present in myPhen
  
  traits = traits[traits %in% names(myPhen)]
  if (is.data.frame(myPhen2)) traits = traits[traits %in% names(myPhen2)]
  
  subtitle = 1

  message('Creating the population partition matrices and subsetting genotypes with phenotypic and genotypic data available.','\n')
  
  for (trait in traits){
   
    # Get the list of lines that have all variables phenotyped
    all_phen_list = rownames(myPhen)[! is.na(myPhen[,trait])]
    
    # Get the list of samples with all phenotyped variables and genotypic data
    all_phen_gen = intersect(as.character(myGen), all_phen_list)
    
    namesList[[trait]]$all_phen_gen = sort(all_phen_gen)
    
    # If there is a 2nd dataset, get the list of lines that will be analyzed from both datasets
    if (is.data.frame(myPhen2)){
      
      
      all_phen_list2 = rownames(myPhen2)[! is.na(myPhen2[,trait])]
      
      all_phen_gen2 = intersect(as.character(myGen), all_phen_list2)
      
      namesList[[trait]]$all_phen_gen2 = all_phen_gen2 # all_phen_gen2 contains all lines with geno and pheno data in the Validation population (including those that are not present in the Training population)
      
      # Get the list of samples with pheno & geno data from both datasets
      
      all_phen_gen = intersect(all_phen_gen, all_phen_gen2)

      namesList[[trait]]$all_phen_gen = sort(all_phen_gen)

      if (! predNewLines){

        # all_phen_gen2 still contains lines with geno and phen data for myPhen2, but namesList[[trait]]$all_phen_gen2 will not.
        namesList[[trait]]$all_phen_gen2 = all_phen_gen
      }
  
      # Get the size for training and validation populations
      
      namesList[[trait]]$len = c('apg1' = length(all_phen_gen), # How many lines have geno and phen data for both datasets
                                 'apg2' = length(all_phen_gen2) - length(all_phen_gen)) # How many lines have geno and phen data and are unique for myPhen2
      
      if(subtitle) message('  ├─\tTrait\t=\tTP-length\tVP-length\n') ; subtitle = 0
      
      message(' ├─ ', trait,'\t=\t',
              round(length(all_phen_gen) * (1 - validation)),'\t',
              round(length(all_phen_gen) * validation) + length(namesList[[trait]]$all_phen_gen2) - length(all_phen_gen))
      
    } else {
      
      # Get the total size of population
      
      if(subtitle) message('  ├─\tTrait\t=\tLines\n') ; subtitle = 0
      
      message(' ├─ ', trait,'\t=\t',length(all_phen_gen))
      namesList[[trait]]$len = c('apg1' = length(all_phen_gen))
      
    }
    
    namesList[[trait]]$combinat = matrix(0, nrow=length(all_phen_gen), ncol=rand_pars, dimnames=list(row=all_phen_gen))
    
    for( i in 1:rand_pars ){
      # 0 = train ; 1 = test
      namesList[[trait]]$combinat[ sample( 1:length(all_phen_gen), length(all_phen_gen) * validation ) , i] = 1
    }
     
  }
  
  return(namesList)

}
