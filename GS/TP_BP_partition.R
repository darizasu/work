

TP_BP_partition = function (phen,X_line,traits,phen2,rand_pars,validation_prop,pnl){
  
  # Object . A function that returns the genotype IDs that have phenotypic and 
  #          genotypic information for each requested trait. Based on this information,
  #          it produces the random assignment matrices for training and validation populations.
  # Input . phen . A data.frame with phenotypic data. First column contains genotype IDs. 
  #         This will be used as Training and Validation, depending on 'phen2'
  # Input . X_line . A character vector with the genotype IDs that are present in the genotypic matrix.
  # Input . traits . A character vector with the traits to be analyzed. These strings are column names in 'phen' and 'phen2'
  # Input . phen2 . A second data.frame (optional) with phenotypic data. This will be used for Validation only
  # Input . rand_pars . Number of random partition populations to be used for prediction ability assessment.
  # Input . validation_prop . A numeric value in (0,1) defining the proportion of the population to be used for validation
  # Input . pnl . Logical. Inlcude lines that were not present in the Training population in the prediction ? Used only when 'phen2' is provided
  # Output . A list with as many items as valid traits there are. 
  #          Each trait is a sublist containing a vector with genotype IDs (lpg), 
  # Authors: dariza and jdelahoz
  #   Last modified: Oct 18, 2019
  
  namesList = list()
  
  # Avoid those traits that are not present in phen
  
  traits = traits[traits %in% names(phen)]
  if (is.data.frame(phen2)) traits = traits[traits %in% names(phen2)]
  
  for (trait in traits){
   
    # Get the list of lines that have all variables phenotyped
    lp = rownames(phen)[! is.na(phen[,trait])]
    
    # Get the list of samples with all phenotyped variables and genotypic data
    lpg = intersect(as.character(X_line), lp)
    
    namesList[[trait]]$lpg = lpg
    
    l2bu = lpg
    
    # If there is a 2nd dataset, get the list of lines that will be analyzed from both datasets
    if (is.data.frame(phen2)){
      
      lp_vp  = rownames(phen2)[! is.na(phen2[,trait])]
      
      lpg_vp = intersect(as.character(X_line), lp_vp)

      if (pnl){

        lpg_vp = setdiff(lpg_vp, lpg)
        l2bu      = union(l2bu, lpg_vp)
        rand_pars = 1

      } else {

        lpg_vp = intersect(lpg, lpg_vp)
        lpg    = lpg_vp
        lpg_vp = 1:(length(lpg_vp) * validation_prop)
        l2bu   = lpg
      }

      namesList[[trait]]$lpg    = sort(lpg)
      namesList[[trait]]$lpg_vp = sort(lpg_vp)
    }

    namesList[[trait]]$l2bu = l2bu # Lines to be used

    namesList[[trait]]$pop_split = matrix(0, nrow=length(l2bu), ncol=rand_pars,
                                         dimnames=list(row=l2bu))
    
    for( i in 1:rand_pars ){
      # 0 = train ; 1 = test
      if (pnl){

        namesList[[trait]]$pop_split [ lpg_vp , i] = 1

      } else {

        namesList[[trait]]$pop_split [ sample( 1:length(lpg), length(lpg) * validation_prop ) , i ] = 1
      } 
    }
  }
  
  return(namesList)
}
