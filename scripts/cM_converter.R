cM_converter <-  
  function(pred, chrom_col = 1, pos_col = 2, what = 'phys2gen',
           map = 'https://raw.githubusercontent.com/darizasu/work/master/scripts/test_data/MGC_marker_placements_for_genetic_map.map'){
  
  #   Goal      : Convert physical to genetic positions (or genetic to physical positions)
  #               fitting a smooth spline with a physical-genetic map.
  #   Input     : pred . A data.frame containing the chromosome and position 
  #               (physical or genetic accordingly) to be converted.
  #   Input     : chrom_col . Integer indicating the column number containing the chromosome to be converted from 'pred'.
  #   Input     : pos_col . Integer indicating the column number containing the positions to be converted from 'pred'.
  #   Input     : what . Character. Use 'phys2gen' if 'pred' contains physical positions 
  #               to be converted to genetic positions.
  #               Use 'gen2phys' if 'pred' contains genetic positions to be converted to physical positions.
  #   Input     : map . An original tab-separated file containing physical and genetic reference
  #               coordinates to be used to fit the smooth spline. This file must have a header. 
  #               Columns are: chromosome - empty_col - genetic_pos - physical_pos
  #   Output    : Original 'pred' data.frame with the column 'PRED' for predicted genetic positions.
  #   Authors   : japaricio and darizasu
  #    Last update: June 7th, 2019
  
  map = read.table(map, col.names=c('CHR','na','GEN','PHY'))
  
  
  colnames(pred)[c(chrom_col,pos_col)] = c('CHR','POS')
  
  for (chr in levels(map$CHR)){
    
    # pred[pred$CHR == chr,'POS'] = sort(pred[pred$CHR == chr,'POS'])
    
    sbMap <- map[map$CHR == chr,]
    sbPred <- pred[pred$CHR == chr,]
    
    fit <- smooth.spline(sbMap$PHY, sbMap$GEN)
    
    if (what == 'phys2gen'){
      
      sbPred <- sort(predict(fit, sbPred$POS)$y)
      
      if (min(sbPred) < 0){
        
        pred[pred$CHR == chr,'PRED'] <- sbPred + abs(min(sbPred))
        
      } else {
        
        pred[pred$CHR == chr,'PRED'] <- sbPred - min(sbPred)
      }
      
    } else if (what == 'gen2phys'){
      
      sbPred <- approx(x = fit$y, y = fit$x, xout = sbPred$POS)$y
      
      pred[pred$CHR == chr,'PRED'] <- round(sbPred)
      
    } else stop("The argument 'what' must be either 'phys2gen' or 'gen2phys'" )
    
  }
  
  return(pred)
}

# pred = cM_converter(pred = 'MGC_marker_placements_for_genetic_map.map', chrom_col = 1, pos_col = 4)
# write.table(pred[,c(1,2,5,4)], 'MGC_marker_placements_for_genetic_map_predicted.map', col.names=F, quote=F, append=F, row.names=F, sep='\t')

