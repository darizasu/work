clusterPairs <- function(data){
  
  
  # Object : Take pairs of identical samples and find their corresponding clusters
  # Input  : data . A data.frame of two columns and one row for each pair of identical samples
  # Output : A matrix with one cluster per row and its corresponding members as columns
  # Authors: darizasu
  #  Last update: December 8, 2018
  
  require(data.table)
  
  # Remove duplicated pairs of samples
  cat('\nRemoving duplicated pairs in the table ...\n')
  data <- t(apply(X = data, MARGIN = 1, sort))
  data <- unique(data)
  
  # In case there is any pair with the same genotype, place it in the bottom of the matrix
  if(any(data[,1] == data[,2])){
    
    orphans <- data[data[,1] == data[,2],]
    data <- data[data[,1] != data[,2],]
    data <- rbind(data, orphans)
  }
  
  # Sort the matrix by 1st and then by 2nd column
  data <- data[do.call(order, as.data.frame(data)),]
  
  # Put all genotypes in a vector to display at the end
  nGenIDs <- length(unique(as.vector(data)))
  
  # Remove any remaining rownames in the matrix
  dimnames(data) = NULL
  
  # This list will contain one element per cluster
  myList = list()
  
  # Set the progress bar
  pb <- txtProgressBar(min = 1, max = nrow(data), initial = 0, style = 3)
  
  # Place the first row of the matrix in the list
  if(length(myList) == 0) myList[[1]] <- unname(data[1,])
  
  cat('Finding clusters for ', nGenIDs, ' samples:\n')
  for(i in 2:nrow(data)){
    
    # For each row, 'A' and 'B' are logical vectors indicating which element(s)
    # of the list contains the first and second elements of the pair of samples in the current row
    # Use 'unname' to drop any names in the vector
    A <- unname(data[i,1]) %chin% myList[[1]]
    B <- unname(data[i,2]) %chin% myList[[1]]
    
    if (length(myList) > 1){
      
      for (j in 2:length(myList)){
        
        A <- c(A, unname(data[i,1]) %chin% myList[[j]])
        B <- c(B, unname(data[i,2]) %chin% myList[[j]])
      }
    }
    
    # Merge positions of 'A' and 'B' into 'C'
    C <- A | B
    
    # In case any of the genotypes in the current row is already present in any cluster, 
    # merge the elements of the cluster with the elements of the row
    if (any(C)){
      
      # If there's only one cluster containing any of the genotypes in the row, 
      # merge the elements of the cluster with the elements of the row and continue.
      # If there are more than 1 clusters containing any of the genotypes in the row,
      # merge those clusters into a single one and then 
      # merge the elements of that new cluster with the elements of the row and continue.
      if(sum(C) == 1){
        
        myList[[which(C)]] <- union(myList[[which(C)]], unname(data[i,]))
        
      } else {
        
        myList[[which(C)[1]]] <- unique(unlist(myList[which(C)]))
        myList[[which(C)[1]]] <- union(myList[[which(C)[1]]], unname(data[i,]))
        myList[which(C)[-1]] <- NULL
      }
      
    # If none of the genotypes is present in any cluster, create a new cluster
    } else {
      
      myList[[length(myList) + 1]] <- unname(data[i,])
    } 
    
    # Notify the current progress
    setTxtProgressBar(pb = pb, value = i)
  }
  
  # Convert myList to a matrix
  myClusters <- t(sapply(myList, "[", seq(max(lengths(myList)))))
  
  # Print a brief summary and close the progress bar
  cat('\n\n', length(myList),'\t\tclusters were identified for ', nGenIDs, ' samples.\n', sep = '')
  
  uniqueClusters <- sum( apply(X = myClusters, MARGIN = 1, function(x) length( unique( na.omit(x) ))) == 1)
  cat(uniqueClusters,'\t\tclusters had only 1 sample.\n', sep = '')
  
  close(pb)
  
  # Return the final matrix with clusters
  return(myClusters)
}
