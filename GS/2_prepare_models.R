 # y = vector of values for the observed trait
 # X = X matrix (samples x markers (-1, 0, 1))
 # pop.split = random permutations of samples
 # model= c( "BayesA", "BayesB", "BayesC", "BayesRR", "BLasso", "BLassof", "FIXED", "RKHS", "GBLUP" )#,"rrBLUP")

setwd("/bioinfo1/projects/bean/VEF/genomic_selection/scripts")

runBGLR <- function(y, trait, X, pop.split, yBP, model, saveAt = paste(outDir,'/',sep=''),
                    Gmatrix = '../geno/VEF_noMeso_annotated_repMasked_q40_s_fi_maf05_oh06_i210_imputed_kinship.txt')
{
	       if (model=="BayesA"){
		ETA = list(list(model="BayesA", X=X))
	} else if (model=="BayesB") {
		ETA = list(list(model="BayesB", X=X))
	} else if (model=="BayesC") {
		ETA = list(list(model="BayesC", X=X))
	} else if (model=="BayesRR") {
		ETA = list(list(model="BRR",    X=X))
	} else if (model=="BLasso") {
		ETA = list(list(model="BL",     X=X))
	} else if (model=="FIXED") {
		ETA = list(list(model="FIXED",  X=X))
	} else if (model=="BLassof") {
		# Z = scale(X)
		# G = tcrossprod(Z) / ncol(Z)
	  G = read.delim(Gmatrix, row.names = 1)
	  G = as.matrix(G[all_phen_gen,all_phen_gen])
		L = svd(G)
		Lm = L$u %*% diag(L$d)^(1/2)
		ETA = list(list(model="BL",     X=Lm))
	} else if (model=="RKHS") {
		D = as.matrix( dist( X, method="euclidean")) ^2
		D = D / mean(D)
		h = 0.5
		K = exp(-h * D)
		ETA = list(list(model="RKHS",   K=K))
	} else if (model=="GBLUP") {
		# Z = scale(X)
		# G = tcrossprod(Z) / ncol(Z)
	  G = read.delim(Gmatrix, row.names = 1)
	  G = as.matrix(G[all_phen_gen,all_phen_gen])
		ETA = list(list(model="RKHS",   K=G))
	} else {
		cat("ERROR: model",model,"is not available. Continuing...\n")
		return(list( result=NA, cor=NA))
	}
  
  # Phenotypes from the Validation population are set to NA
  
  yNA = y[,trait]
  iNA = which(pop.split == 1)
  yNA[iNA] = NA
  iNA = is.na(yNA) # Those lines that were not phenotyped in the 1st dataset will be predicted as well.
  
  fm = BGLR( y=yNA, ETA=ETA, 
             nIter=1000, burnIn=200, thin=10, 
             verbose=FALSE, saveAt=saveAt)
  
  TP = fm$yHat[iNA] # Predicted phenotypes from validation population
    
  if (missing(yBP)){
    
    BP = y[iNA,trait]
    
  } else {
    
    BP = yBP[iNA,trait]
    
  }
	
  return( list(result = fm,
               cor    = cor(TP, BP, use="complete.obs", method="pearson")  ) )
  
}
