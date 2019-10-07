intertek2rQTL <- 
  function(dta, P1, P2, markerID, P1c, P2c){
  
    #   Goal      : Convert a genotypic matrix from the Intertek format to the R/qtl2 format.
    #   Input     : dta . A data.frame containing the genotypic matrix in the intertek format.
    #               In this format, the first column is the marker names. The first row is the individual IDs.
    #               There are one column for each founder line with genotype calls in the form 'A/A', 'C/C', etc.
    #   Input     : P1 . Character indicating the column name of the first founder.
    #   Input     : P2 . Character indicating the column name of the second founder.
    #   Input     : markerID . Character indicating the column name of the marker IDs (The first column).
    #   Input     : P1c . Character indicating the letter to use for the genotype calls derived from the founder P1.
    #               This is usually the first letter of the founder P1. Make sure it is different from P2c.
    #   Input     : P2c . Character indicating the letter to use for the genotype calls derived from the founder P2.
    #               This is usually the first letter of the founder P2. Make sure it is different from P1c.
    #   Output    : A data.frame with the genotypic matrix in the R/qtl2 format.
    #               The first column is the individual IDs and the first row is the marker names.
    #   Authors   : darizasu
    #    Last update: October 7th, 2019
  
  # Make sure the table is a data.frame, not a tibble
  dta <- as.data.frame(dta)
  
  # Remove markers with no ID, or markers with missing calls in any founder
  noID <- is.na(dta[,markerID]) | is.na(dta[,P1]) | is.na(dta[,P2])
  
  if (sum(noID)){
    
    dta <- dta[ ! noID , ]
    message('\n', sum(noID), ' marker(s) was/were removed from the dataset because of missing calls for any of the founders or no markerID.\n')
  }
  
  # Find duplicated rows and remove them
  dups <- duplicated(dta)
  
  if (sum(dups)){
    
    message('The marker(s) ID(s):\n', paste(dta[dups,markerID], sep = '\n'),'\nis/are duplicated in the original table.\n')
    dta <- dta[ !dups, ]
  }
  
  # Find duplicated marker IDs. Stop if any.
  dups <- duplicated(dta[,markerID])
  
  if (sum(dups)){
    
    stop('The marker(s) ID(s):\n\n', dta[dups,markerID], '\n\nis/are duplicated in the original table, but they have different genotype calls in the population.\n')
  }
  
  # Replace 'A/A' to 'A' in the founders' alleles
  dta[,P1] <- gsub(pattern = '([[:alpha:]])\\1+', replacement = '\\1', 
                   gsub(pattern = '/', replacement = '', dta[,P1]))
  
  dta[,P2] <- gsub(pattern = '([[:alpha:]])\\1+', replacement = '\\1', 
                   gsub(pattern = '/', replacement = '', dta[,P2]))
  
  # Split the table into a list by markers, remove unused levels of factors
  dta <- split(x = dta, f = dta[,markerID], drop = T)
  dta <- lapply(dta, droplevels)
  
  # For each marker ...
  dta <- mapply(function(marker){
    
    # Make a vector of values to be replaced.
    # The names of the vector are replacement values.
    # This works as a key:value pair
    mp <- setNames(c(marker[,markerID],     P1c,        P2c),
                   c(marker[,markerID], marker[P1], marker[P2]))

    # Get a vector of unique values for this marker to check how many different genotypes are in it.
    uniq <- unique( unlist(marker) )
    
    # Find heterozygous calls not present in the vector of the values to be replaced (mp). Do not include NAs
    het <- uniq[ ! uniq %in% names(mp)]
    het <- het[!is.na(het)]
    
    # In case this marker has heterozygous calls
    if (length(het)){
      
      # Create a vector with the missing replacement values of 'mp'
      het <- setNames(sapply(strsplit(x = het, split = '/'), 
                             function(x){ paste0(mp[x], collapse = '') }),
                      het)
      
      # Add the missing value to 'mp'
      mp <- c(mp, het)
    }

    # Change Intertek format to rQTL format for this marker
    marker[] <- mp[unlist(marker)]
    
    # Return this marker as an element of the list
    return(marker)
    
  }, dta, SIMPLIFY = F, USE.NAMES = T)
  
  # Pile up the elements of the list to get a data.frame
  dta <- do.call("rbind", dta)
  
  # Remove rownames and transpose the dataframe for rQTL format
  rownames(dta) <- NULL
  dta <- setNames( data.frame(t(dta[,-1])), dta[,1] )
  
  return(dta)
  
}


# mp <- setNames(c(dta$`BARC-PV-0000055`[,markerID], P1c, P2c),
#                c(dta$`BARC-PV-0000055`[,markerID], dta$`BARC-PV-0000055`[P1], dta$`BARC-PV-0000055`[P2]))
# 
# uniq <- unique(unlist(dta$`BARC-PV-0000055`))
# 
# if (length(uniq) > 4 & length(uniq) <= 6){
#   
#   het <- uniq[ ! uniq %in% names(mp)]
#   het <- het[!is.na(het)]
#   
#   het <- setNames(sapply(strsplit(x = het, split = '/'),
#                          function(x){ paste0(mp[x], collapse = '') }),
#                   het)
#   
#   mp <- c(mp, het)
#   
# }
# 
# dta$`BARC-PV-0000055`[] <- mp[unlist(dta$`BARC-PV-0000055`)]
# 
