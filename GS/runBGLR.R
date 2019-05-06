
runBGLR <- function(y, trait, X, pop.split, yBP, model, G, saveAt = paste(outDir,'/',sep=''), myNames){
  
  # Object . A function that handles input objects to run BGLR
  # Input . y . Vector of values for the observed trait
  # Input . trait . 
  # Input . X . Genotypic matrix (samples x markers (-1, 0, 1))
  # Input . pop.split . Random permutations of samples
  # Input . yBP .
  # Input . model . c( "BayesA", "BayesB", "BayesC", "BayesRR", "BLasso", "BLassof", "FIXED", "RKHS", "GBLUP")
  # Input . G .
  # Input . saveAt .
  # Input . myNames .
  # Output . A list with fitted model and the correlation between observed and predicted values.
  # Authors: jfdelahoz and darizasu
  #   Last modified: November 2, 2018
  
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
    
    L = svd(G)
    Lm = L$u %*% diag(L$d)^(1/2)
    ETA = list(list(model="BL",     X=Lm))
    
  } else if (model=="RKHS") {
    
    D = as.matrix( dist( X[myNames,], method="euclidean")) ^2
    D = D / mean(D)
    h = 0.5
    K = exp(-h * D)
    ETA = list(list(model="RKHS",   K=K))
    
  } else if (model=="GBLUP") {
    
    ETA = list(list(model="RKHS",   K=G))
    
  } else {
    
    stop("Model",model,"is not available.")
  }
  
  # Phenotypes from the Validation population are set to NA
  
  yNA = y[myNames,trait]
  iNA = which(pop.split == 1)
  yNA[iNA] = NA
  iNA = is.na(yNA) # Those lines that were not present in the Training population will be predicted as well.
  
  fm = BGLR( y=yNA, ETA=ETA, 
             nIter=10000, burnIn=1000, thin=5, 
             verbose=FALSE, saveAt=saveAt)
  
  TP = fm$yHat[iNA] # Predicted phenotypes from validation population
  
  if (missing(yBP)){
    
    BP = y[iNA,trait]
    
  } else {
    
    yBP = yBP[myNames,]
    BP = yBP[iNA,trait]
    
  }
  
  return( list(result = fm,
               cor    = cor(TP, BP, use="pairwise.complete.obs", method="pearson")  ) )
}
