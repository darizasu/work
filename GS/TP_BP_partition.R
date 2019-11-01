

TP_BP_partition = function (myPhen,myGen,traits,myPhen2,rand_pars,validation,pnl){
  
  # Object . A function that returns the genotype IDs that have phenotypic and 
  #          genotypic information for each requested trait. Based on this information,
  #          it produces the random assignment matrices for training and validation populations.
  # Input . myPhen . A data.frame with phenotypic data. First column contains genotype IDs. 
  #         This will be used as Training and Validation, depending on 'phen2'
  # Input . myGen . A character vector with the genotype IDs that are present in the genotypic matrix.
  # Input . traits . A character vector with the traits to be analyzed. These strings are column names in 'phen' and 'phen2'
  # Input . myPhen2 . A second data.frame (optional) with phenotypic data. This will be used for Validation only
  # Input . rand_pars . Number of random partition populations to be used for prediction ability assessment.
  # Input . validation . A numeric value in (0,1) defining the proportion of the population to be used for validation
  # Input . pnl . Logical. Inlcude lines that were not present in the Training population in the prediction ? Used only when 'phen2' is provided
  # Output . A list with as many items as valid traits there are. 
  #          Each trait is a sublist containing a vector with genotype IDs (linePG), 
  # Authors: dariza and jdelahoz
  #   Last modified: Oct 18, 2019
  
  namesList = list()
  
  # Avoid those traits that are not present in myPhen
  
  traits = traits[traits %in% names(myPhen)]
  if (is.data.frame(myPhen2)) traits = traits[traits %in% names(myPhen2)]
  
  for (trait in traits){
   
    # Get the list of lines that have all variables phenotyped
    lineP = rownames(myPhen)[! is.na(myPhen[,trait])]
    
    # Get the list of samples with all phenotyped variables and genotypic data
    linePG = intersect(as.character(myGen), lineP)
    
    namesList[[trait]]$linePG = linePG
    
    l2bu = linePG
    
    # If there is a 2nd dataset, get the list of lines that will be analyzed from both datasets
    if (is.data.frame(myPhen2)){
      
      lineP_vp  = rownames(myPhen2)[! is.na(myPhen2[,trait])]
      
      linePG_vp = intersect(as.character(myGen), lineP_vp)

      if (pnl){

        linePG_vp = setdiff(linePG_vp, linePG)
        l2bu      = union(l2bu, linePG_vp)
        rand_pars = 1

      } else {

        linePG_vp = intersect(linePG, linePG_vp)
        linePG    = linePG_vp
        linePG_vp = 1:(length(linePG_vp) * validation)
        l2bu      = linePG
      }

      namesList[[trait]]$linePG    = sort(linePG)
      namesList[[trait]]$linePG_vp = sort(linePG_vp)
    }

    namesList[[trait]]$l2bu = l2bu # Lines to be used

    namesList[[trait]]$combinat = matrix(0, nrow=length(l2bu), ncol=rand_pars,
                                         dimnames=list(row=l2bu))
    
    for( i in 1:rand_pars ){
      # 0 = train ; 1 = test
      if (pnl){

        namesList[[trait]]$combinat[ linePG_vp , i] = 1

      } else {

        namesList[[trait]]$combinat[ sample( 1:length(linePG), length(linePG) * validation ) , i ] = 1
      } 
    }
  }
  
  return(namesList)
}
