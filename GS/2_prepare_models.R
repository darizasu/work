 # y = vector of values for the observed trait
 # X = X matrix (samples x markers (-1, 0, 1))
 # pop.split = random permutations of samples
 # model= c( "BayesA", "BayesB", "BayesC", "BayesRR", "BLasso", "BLassof", "FIXED", "RKHS", "GBLUP" )#,"rrBLUP")

runBGLR <- function(y, trait, X, pop.split, yBP, model, G, saveAt = paste(outDir,'/',sep=''), myNames)
{
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
		# Z = scale(X[myNames,])
		# G = tcrossprod(Z) / ncol(Z)
	  G = read.delim(G, row.names = 1, header = F)
	  colnames(G) = rownames(G)
	  G = as.matrix(G[myNames,myNames])
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
		# Z = scale(X[myNames,])
		# G = tcrossprod(Z) / ncol(Z)
	  G = read.delim(G, row.names = 1, header = F)
	  colnames(G) = rownames(G)
	  G = as.matrix(G[myNames,myNames])
		ETA = list(list(model="RKHS",   K=G))
	} else {
		cat("ERROR: model",model,"is not available. Continuing...\n")
		return(list( result=NA, cor=NA))
	}
  
  # Phenotypes from the Validation population are set to NA
  
  yNA = y[myNames,trait]
  iNA = which(pop.split == 1)
  yNA[iNA] = NA
  iNA = is.na(yNA) # Those lines that were not phenotyped in the 1st dataset will be predicted as well.
  
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
