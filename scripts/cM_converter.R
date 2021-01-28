cM_converter <- 
  function(markers, chrom_col = 1, pos_col = 2, what = 'phys2gen',
           map = 'https://raw.githubusercontent.com/darizasu/work/master/scripts/test_data/MGC_marker_placements_for_genetic_map.map'){
    
    #   Goal      : Convert physical to genetic positions (or genetic to physical positions)
    #               fitting a smooth spline with a physical-genetic map.
    #   Input     : markers . A data.frame containing the chromosome and position 
    #               (physical or genetic accordingly) to be converted.
    #   Input     : chrom_col . Integer indicating the column number containing the chromosome to be converted from 'markers'.
    #   Input     : pos_col . Integer indicating the column number containing the positions to be converted from 'markers'.
    #   Input     : what . Character. Use 'phys2gen' if 'markers' contains physical positions 
    #               to be converted to genetic positions.
    #               Use 'gen2phys' if 'markers' contains genetic positions to be converted to physical positions.
    #   Input     : map . An original tab-separated file containing physical and genetic reference
    #               coordinates to be used to fit the smooth spline. This file must have a header. 
    #               Columns are: chromosome - empty_col - genetic_pos - physical_pos
    #   Output    : Original 'markers' data.frame with the column 'PRED' for predicted genetic positions.
    #   Authors   : japaricio and darizasu
    #    Last update: January 27th, 2021
    
    map <- read.table(map, col.names=c('Chromosome','na','GEN','PHY'))
    
    colnames(markers)[c(chrom_col,pos_col)] <- c('Chromosome','Position')
    
    for (chr in unique(map$Chromosome)){
      
      sbMap <- map[map$Chromosome == chr,]
      sbPred <- markers[markers$Chromosome == chr,]
      
      fit <- smooth.spline(sbMap$PHY, sbMap$GEN)
      
      if (what == 'phys2gen'){
        
        sbPred <- predict(fit, sbPred$Position)$y
        
        markers[markers$Chromosome == chr,'Predicted_cM'] <- sbPred
        colnames(markers)[pos_col] <- 'Position_bp'
        
      } else if (what == 'gen2phys'){
        
        sbPred <- approx(x = fit$y, y = fit$x, xout = sbPred$Position)$y
        
        markers[markers$Chromosome == chr,'Predicted_bp'] <- round(sbPred)
        
        colnames(markers)[pos_col] <- 'Position_cM'
        
      } else stop("The argument 'what' must be either 'phys2gen' or 'gen2phys'" )
      
    }
    
    return(markers)
  }
