
BGLRwrap <- function(phen, trait, X, pop.split, phen2, model, G, saveAt, myNames){
  
  # Object . A function that creates input objects to run different priors in BGLR
  # Input . phen . A data.frame with phenotypic data. First column contains genotype IDs. This will be used as Training and Validation, depending on 'phen2'
  # Input . trait . A string with the trait to be analyzed. This string is a column name of 'phen'
  # Input . X . Genotypic matrix (samples x markers (-1, 0, 1))
  # Input . pop.split . A numeric vector with random permutations of samples. 0 for training, 1 for validation
  # Input . phen2 . A second data.frame (optional) with phenotypic data. This will be used for Validation only
  # Input . model . c( "BayesA", "BayesB", "BayesC", "BayesRR", "BLasso", "BLassof", "FIXED", "RKHS", "GBLUP")
  # Input . G . A kinship matrix
  # Input . saveAt . A string for the location of the output files.
  # Input . myNames . A character vector with genotype IDs to be used from 'phen' and 'phen2'
  # Output . A list with the fitted model, the correlation between observed and predicted values and a correlation in a data.frame object.
  # Authors: jfdelahoz and darizasu
  #   Last modified: October 18, 2019
  
  if (model=="BayesA"){
    
    ETA = list(list(model="BayesA", X=X[myNames,]))
    
  } else if (model=="BayesB") {
    
    ETA = list(list(model="BayesB", X=X[myNames,]))
    
  } else if (model=="BayesC") {
    
    ETA = list(list(model="BayesC", X=X[myNames,]))
    
  } else if (model=="BayesRR") {
    
    ETA = list(list(model="BRR",    X=X[myNames,]))
    
  } else if (model=="BLasso") {
    
    ETA = list(list(model="BL",     X=X[myNames,]))
    
  } else if (model=="FIXED") {
    
    ETA = list(list(model="FIXED",  X=X[myNames,]))
    
  } else if (model=="BLassof") {
    
    L = svd(G[myNames,myNames])
    Lm = L$u %*% diag(L$d)^(1/2)
    ETA = list(list(model="BL",     X=Lm))
    
  } else if (model=="RKHS") {
    
    D = as.matrix( dist( X[myNames,], method="euclidean")) ^2
    D = D / mean(D)
    h = 0.5
    K = exp(-h * D)
    ETA = list(list(model="RKHS",   K=K))
    
  } else if (model=="GBLUP") {
    
    ETA = list(list(model="RKHS",   K=G[myNames,myNames]))
    
  } else {
    
    stop("Model",model,"is not available.")
  }
  
  # Phenotypes from the Validation population are set to NA
  
  yNA = phen[myNames,trait]
  iNA = which(pop.split == 1)
  yNA[iNA] = NA
  
  fm = BGLR( y=yNA, ETA=ETA, 
             nIter=10000, burnIn=1000, thin=5, 
             verbose=FALSE, saveAt=saveAt)
  
  Pred = fm$yHat[iNA] # Predicted phenotypes from validation population
  
  if (missing(phen2)){
    
    Obs = phen[myNames,trait]
    Obs = Obs[iNA]
    
  } else {
    
    Obs = phen2[myNames,trait]
    Obs = Obs[iNA]
    
  }
  
  return( list(result = fm,
               cor    = cor(Pred, Obs, use="pairwise.complete.obs", method="pearson"),
               table  = data.frame(model = model, trait = trait, pred = Pred, obs = Obs)  ) )
}
