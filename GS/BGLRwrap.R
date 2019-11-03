
BGLRwrap <- function(phen, trait, X, pop_split, phen2, prior, G, out_dir, l2bu, ETA){

  # Object . A function that creates input objects to run different priors in BGLR
  # Input . phen . A data.frame with phenotypic data. First column contains genotype IDs. This will be used as Training and Validation, depending on 'phen2'
  # Input . trait . A string with the trait to be analyzed. This string is a column name of 'phen'
  # Input . X . Genotypic matrix (samples x markers (-1, 0, 1))
  # Input . pop_split . A numeric vector with random permutations of samples. 0 for training, 1 for validation
  # Input . phen2 . A second data.frame (optional) with phenotypic data. This will be used for Validation only
  # Input . prior . c( "BayesA", "BayesB", "BayesC", "BayesRR", "BLasso", "BLassof", "FIXED", "RKHS", "GBLUP")
  # Input . G . A kinship matrix
  # Input . out_dir . A string for the location of the output files.
  # Input . l2bu . A character vector with genotype IDs to be used from 'phen' and 'phen2'
  # Output . A list with the fitted prior, the correlation between observed and predicted values and a correlation in a data.frame object.
  # Authors: jfdelahoz and darizasu
  #   Last modified: October 18, 2019

  if (prior=="BayesA"){

    ETA = list(list(model="BayesA", X=X[l2bu,]))

  } else if (prior=="BayesB") {

    ETA = list(list(model="BayesB", X=X[l2bu,]))

  } else if (prior=="BayesC") {

    ETA = list(list(model="BayesC", X=X[l2bu,]))

  } else if (prior=="BayesRR") {

    ETA = list(list(model="BRR",    X=X[l2bu,]))

  } else if (prior=="BLasso") {

    ETA = list(list(model="BL",     X=X[l2bu,]))

  } else if (prior=="FIXED") {

    ETA = list(list(model="FIXED",  X=X[l2bu,]))

  } else if (prior=="BLassof") {

    ETA = ETA

  } else if (prior=="RKHS") {

    ETA = ETA

  } else if (prior=="GBLUP") {

    ETA = list(list(model="RKHS",   K=G[l2bu,l2bu]))

  } else {

    stop("Model ", prior, " is not available.")
  }

  # Phenotypes from the Validation population are set to NA

  yNA = phen[l2bu,trait]
  iNA = which(pop_split == 1)
  yNA[iNA] = NA

  fm = BGLR( y = yNA, ETA = ETA, 
             nIter = 10000, burnIn = 1000, thin = 5, 
             verbose = FALSE, saveAt = out_dir)

  Pred = fm$yHat[iNA] # Predicted phenotypes from validation population

  if (missing(phen2)){

    Obs = phen[l2bu,trait]
    Obs = Obs[iNA]

  } else {

    Obs = phen2[l2bu,trait]
    Obs = Obs[iNA]
  }

  return( list(result = fm,
               cor    = cor(Pred, Obs, use = "pairwise.complete.obs", method = "pearson"),
               table  = data.frame(model = prior, trait = trait, pred = Pred, obs = Obs)  ) )
}
