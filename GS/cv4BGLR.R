cv4BGLR <- 
  function(model, phen, traits, out_dir = getwd(), save_table = NA, 
           geno, samp, Gmatrix = NA, validation_perc = 30, rand_pars = 10,
           names_list = NA, phen2 = NA, common = TRUE, pnl = FALSE, rand_SNPs = NA){

    # Load BGLR
    if (!requireNamespace("BGLR", quietly = TRUE)) install.packages("BGLR")
    library(BGLR, quietly = T)

    # If pnl, avoid performing cross validation
    rand_pars = ifelse(pnl, 1, rand_pars)

    # Check the training-validation is between [0,100]%, and make it decimal
    if (any(validation_perc >= 100 | validation_perc < 0)){

      stop('The validation population partition should be a value between 0 - 100.')

    } else {

      validation_prop = validation_perc / 100
    }

    # Read geno and pheno data
    X_line = read.delim(samp, header = F, stringsAsFactors = F)[,1]
    X      = read.table(geno, row.names = as.character(X_line), header = F)
    phen   = read.delim(phen, row.names = 1, check.names = F)

    phen2 = if (!is.na(phen2)) read.delim(phen2, row.names = 1, check.names = F) else NA

    # Import or calculate a kinship matrix
    if (is.na(Gmatrix)){

      warning("Using the function 'A.mat' from the package 'rrBLUP' to calculate the kinship matrix.\n", immediate. = T)
      source('https://raw.githubusercontent.com/cran/rrBLUP/master/R/A.mat.R')
      G = A.mat(X)

    } else {

      message("Using the file ", Gmatrix," as the kinship matrix for the population.\n")
      G = read.delim(Gmatrix, row.names = 1, header = F)
      colnames(G) = rownames(G)
    }

    # Make a subset of SNPs if requested. Write the subsetted SNPs to a file for future validation.
    if (!is.na(rand_SNPs)){

      if (ncol(X) < as.integer(rand_SNPs)){
        stop('The number of SNPs to select is higher than the actual number of SNPs available')
      }

      X = X[,sample(1:ncol(X), as.integer(rand_SNPs))]

      write.table(X, paste0(out_dir,'/Geno_matrix_',rand_SNPs,'SNPs_rrBLUP.in'),
                  sep = ' ', col.names = F, quote = F, row.names = F)

      message(rand_SNPs," SNPs have been selected from the original genotype matrix.",
              " The selected SNPs were saved at ",
              paste0(out_dir,'/Geno_matrix_',rand_SNPs,'SNPs_rrBLUP.in'),'\n')
    }

    # Import complementary functions
    source("https://raw.githubusercontent.com/darizasu/work/master/GS/TP_BP_partition.R")
    source("https://raw.githubusercontent.com/darizasu/work/master/GS/BGLRwrap.R")

    # Generate partitions Training-Validation according to the input data provided.
    if (is.na(names_list)){

      message('Creating the population partition matrices and subsetting genotypes with phenotypic and genotypic data available.\n')

      names_list = TP_BP_partition(phen = phen, X_line = X_line, traits = traits,
                                   phen2 = phen2, rand_pars = rand_pars, validation_prop = validation_prop, pnl = pnl)

      save(names_list, file = paste0(out_dir,'/names_list.RData'))
      message('The population partition matrices and genotype names have been saved at: ',
              paste0(out_dir,'/names_list.RData'),'\n')

      # Import partitions Training-Validation if provided  
    } else {

      nlf = names_list # The string containing the location of the .Rdata file is saved in nlf
      message('Loading population partition matrices and genotype names from:',names_list,'\n')
      load(names_list) # The 'name_list' object is no longer a string with the location of the .Rdata file, instead is the list that was contained in the file.

      for (trait in names(names_list)){

        if (ncol(names_list[[trait]]$pop_split) < rand_pars){
          stop(paste('There are only',ncol(names_list[[trait]]$pop_split),
                     'population partitions available at',nlf,
                     'for the trait',trait,'and you requested',rand_pars,
                     "partitions with '-i' option."))
        }
      }
    }

    # Find traits with the same set of genotypes and assign a common Training-Validation partition.
    if (common){

      rand_pops = lapply(names_list, `[[`, 'l2bu')
      unique_pops = unique(rand_pops)

      pops_n = sapply(rand_pops, function(x) which( sapply(unique_pops, identical, x) ))

      message('\n\nA common population partition matrix will be assigned to traits with exactly the same set of genotypes:\n')

      for (i in unique(pops_n)){

        pop = names(pops_n[pops_n == i])
        first = names_list[[ pop[1] ]]$pop_split
        message(paste0(' |-- ', pop, coll = '\n'))

        for (j in pop) names_list[[ j ]]$pop_split = first
      }
    }

    # Display information about the Training-Validation composition for each trait
    sbt = 1
    message('\n\nGenotype composition for each trait:\n')

    for (trait in traits){

      if (is.data.frame(phen2)){

        if (sbt) message('  --\tTrait\t=\tTP-length\tVP-length\n')
        message(' -- ',trait,'\t=\t', length(names_list[[trait]]$lpg),'\t=\t',
                length(names_list[[trait]]$lpg_vp))
        sbt = 0

      } else {

        if (sbt) message('  --\tTrait\t=\tLines\n')
        message('  --  ',trait,'\t=\t',length(names_list[[trait]]$lpg))
        sbt = 0
      }
    }

    cat('\nmodel\ttrait\trandomPop\tcorr\tfinishedAt\n')

    out_table = data.frame(prior = NA, trait = NA, randomPop = NA, corr = NA, finishedAt = NA)

    priors = setNames(c('BayesA','BayesB','BayesC','BRR',    'BL',    'FIXED'),
                      c('BayesA','BayesB','BayesC','BayesRR','BLasso','FIXED'))

    # Run cross validation
    for (trait in traits){

      set.seed(1234)
      pop_split = names_list[[trait]]$pop_split
      l2bu  = names_list[[trait]]$l2bu

      for (prior in model){

        if(prior == "BLassof"){

          L = svd(G[l2bu,l2bu])
          Lm = L$u %*% diag(L$d) ^ (1/2)
          ETA = list(list(model = "BL", X = Lm))

        } else if(prior == "RKHS"){

          D = as.matrix( dist( X[l2bu,], method="euclidean")) ^ 2
          D = D / mean(D)
          h = 0.5
          K = exp(-h * D)
          ETA = list(list(model = "RKHS", K = K))

        } else if(prior == "GBLUP"){

          ETA = list(list(model = "RKHS",   K = G[l2bu,l2bu]))

          } else ETA = list(list(model = priors[prior], X = X[l2bu,]))

        for (i in 1:rand_pars){

          if (! is.data.frame(phen2)){

            if (trait %in% names(phen)){

              out_cor = BGLRwrap(phen = phen,
                               trait = trait,
                               pop_split = pop_split[l2bu,i],
                               prior = prior,
                               l2bu = l2bu,
                               out_dir = paste0(out_dir,'/'),
                               ETA = ETA)

              cat(prior, trait, paste0('pop',i), out_cor$cor, paste0(Sys.time(),'\n'), sep = '\t')

              out_table = rbind(out_table, data.frame(prior = prior, 
                                                  trait = trait, 
                                                  randomPop = paste0('pop',i), 
                                                  corr = out_cor$cor, 
                                                  finishedAt = Sys.time()))
            }

          } else {

            shared_traits = intersect(names(phen), names(phen2))

            if (trait %in% shared_traits){

              out_cor = BGLRwrap(phen = phen,
                               trait = trait,
                               pop_split = pop_split[l2bu,i],
                               phen2 = phen2,
                               prior = prior,
                               l2bu = l2bu,
                               out_dir = paste0(out_dir,'/'),
                               ETA = ETA)

              cat(prior, trait, paste0('pop',i), out_cor$cor, paste0(Sys.time(),'\n'), sep = '\t')

              out_table = rbind(out_table, data.frame(prior = prior, 
                                                  trait = trait, 
                                                  randomPop = paste0('pop',i), 
                                                  corr = out_cor$cor, 
                                                  finishedAt = Sys.time()))
            }
          }
        }
      }
    }

    if(!is.na(save_table)) write.csv(out_table, save_table, quote = F, row.names = F)
    return(out_table)
  }
