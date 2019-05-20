ggCor <- 
  function(myData, colours = c('#db4437','white','#FF9D00'),
           blackLabs = c(-0.7, 0.7), showSignif = TRUE,
           pBreaks = c(0, .0001, .001, .01, Inf), pLabels = c('***','**','*', 'ns'),
           showDiagonal = FALSE, Diag = NULL, returnTable = FALSE){
    
  #   Goal      : Return a ggplot object to plot a triangular correlation figure between 2 or more variables.
  #               Depends on the packages 'ggplot2' 'psych' and 'reshape'
  #
  #   Input     : myData       = A data.frame with numerical columns for each variable to be compared.
  #   Input     : colours      = A vector of size three with the colors to be used for values -1, 0 and 1.
  #   Input     : blackLabs    = A numeric vector of size two, with min and max correlation coefficient 
  #                              limits to display with black tags. Any value outside this range will be 
  #                              displayed with white tags.
  #   Input     : showSignif   = Logical scalar. Display significance values ?
  #   Input     : pBreaks      = Passed to function 'cut'. Either a numeric vector of two or more unique 
  #                              cut points or a single number (greater than or equal to 2) giving the
  #                              number of intervals into which x is to be cut.
  #   Input     : pLabels      = Passed to function 'cut'. labels for the levels of the resulting category.
  #                              By default, labels are constructed using "(a,b]" interval notation. 
  #                              If pLabels = FALSE, simple integer codes are returned instead of a factor.
  #   Input     : showDiagonal = Logical scalar. Display main diagonal values ?
  #   Input     : Diag         = A named vector of labels to display in the main diagonal. The names are 
  #                              used to place each value in the corresponding coordinates of the diagonal.
  #                              Hence, these names must be the same as the colnames of myData
  #   Input     : returnTable  = Return the table to display instead of a ggplot object
  #
  #   Output    : A ggplot object containing a triangular correlation figure with all numeric variables 
  #               in myData. If returnTable is TRUE, the table used to produce the figure is returned instead.
  #   Authors   : darizasu
  #    Last update: May 18, 2019
    
  
  # Drop non numeric columns in the dataset
  if (sum( !sapply(myData, is.numeric) )){
    
    message('Dropping non-numeric columns in the dataset:\n',
            paste(names( which(!sapply(myData, is.numeric)) ),
                  collapse = '\t'))
    
    myData = myData[,sapply(myData, is.numeric)]
  }
  
  # Calculate corr-coeffs and p values
  cors = psych::corr.test(myData, use = 'pairwise.complete.obs')
  
  # Use the adjusted p values for multiple testing instead of raw coeffs
  cors$p = t(cors$p)
  
  # Keep only the matrices with correlation coefficients and p values
  cors = cors[c(1,4)]
  
  # For each matrix, do ...
  cors = lapply(cors, function(x){
    
    # Keep the upper triangle of the matrix
    x[upper.tri(x)] = NA
    
    # Transpose the matrix to plot the lower triangle
    x  = as.data.frame(t(x))
    
    # Reshape the matrix to tidy format
    x[,'col'] = colnames(x)
    x  = reshape::melt(x, id='col')
    colnames(x) = c('col','row','value')
    
    # Round coefficients
    x$name = round(x$value,2)
    
    # Sort the x axis according to myData column order
    x$col = factor(x$col, levels = colnames(myData))
    
    # Reverse the y axis for a triangle plot from top-left to bottom-right
    x$row = factor(x$row, levels = rev(colnames(myData)))
    
    # Remove NAs
    x = na.omit(x)
    
  })
  
  # Combine both dataframes with p values and corr coefficients
  cors = merge(x = cors$r, y = cors$p, by = c('col','row'))
  
  # Keep x, y, p val and corr-coefficients columns
  cors = cors[,c(1,2,4,5)]
  
  if (showSignif){
    
    # Create a categorical variable for p values as defined by pBreaks
    cors$signi = cut(x = cors$value.y,  right = F,
                     breaks = pBreaks, labels = pLabels)
    
    # Join corr-coeff and p-value to display it as a label for each tile
    cors$label = paste(cors$name.x, cors$sign, sep='\n')
    
  } else {
    
    # The label for each tile is the corr-coeff only
    cors$label = cors$name.x
  }
  
  # If there are user-specified values to display in the diagonal
  if (! is.null(Diag)){
    
    # Check the names in Diag are the same than colnames of myData
    if ( sum(! names(Diag) %in% colnames(myData)) ){
      warning("These elements in 'Diag' do not correspond to column names in 'myData':\n",
              paste(names(Diag)[!names(Diag) %in% colnames(myData)],
                    collapse = '\t'))
    }
    
    # The tiles of the diagonal are gray
    cors[cors$col == cors$row, 'name.x'] = NA
    
    # Get the name of x and y levels
    d = as.character(cors[cors$col == cors$row, 'row'])
    
    # Modify the elements of the diagonal and make sure they are displayed
    cors[cors$col == cors$row, 'label'] = Diag[d]
    showDiagonal = TRUE
  }
  
  # Remove the elements of the main diagonal if you don't want to display
  if (!showDiagonal)  cors = cors[cors$col != cors$row,]
  
  # Show darker tiles with white labels for clarity
  cors$txtCol = ifelse(cors$name.x > blackLabs[1] & 
                         cors$name.x < blackLabs[2], 'black', 'white')
  
  # Do not show tile labels for empty tiles.
  # Make tile labels of the diagonal white
  cors$txtCol[is.na(cors$txtCol)] = 'white'
  
  if (returnTable) return(cors)
  
  require(ggplot2)
  
  p = ggplot(data = cors, aes(x = col, y = row, fill = name.x)) + 
    geom_tile(color = 'gray') + labs(x = NULL, y = NULL) + theme_minimal(base_size = 13) +
    geom_text(aes(x = col, y = row, label = label), color = cors$txtCol, size = 3) +
    scale_fill_gradient2(low = colours[1], mid = colours[2], high = colours[3]) + 
    theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.position = 'none',
          panel.grid.minor.x = element_blank(), panel.grid.major = element_blank())
  
  return(p)
}
