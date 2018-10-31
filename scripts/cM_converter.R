cM_converter <-  
  function(map, pred, chrom_col = 1, pos_col = 2, what = 'phys2gen'){
  
  # Object : Convert physical to genetic positions fitting a smooth spline with a physical-genetic map.
  # Input  : map . An original tab-separated physical-genetic map with a header. 
  #          Columns are chromosome - empty_col - genetic_pos - physical_pos
  # Input  : pred . A data.frame containing the chromosome and physical_pos to be converted to genetic_pos
  # Input  : chrom_col . Integer indicating position for chromosome column
  # Input  : map . Integer indicating position for physical position to be converted
  # Output : Original 'pred' data.frame with the column 'GEN' for predicted genetic positions.
  # Authors: japaricio and darizasu
  #  Last update: August 23, 2018
  
  map = read.table(map, col.names=c('CHR','na','GEN','PHY'))
  
  initial = read.table(file=pred, nrows=100, comment.char='', header=T)
  classes = sapply(initial, class)
  pred = read.table(file=pred, colClasses=classes, header=T)
  
  colnames(pred)[c(chrom_col,pos_col)] = c('CHR','POS')
  
  for (chr in levels(map$CHR)){
    
    sbMap <- map[map$CHR == chr,]
    sbPred <- pred[pred$CHR == chr,]
    
    fit <- smooth.spline(sbMap$PHY, sbMap$GEN)
    
    if (what == 'phys2gen'){
      
      sbPred <- sort(predict(fit, sbPred$POS)$y)
      
      if (min(sbPred) < 0){
        
        pred[pred$CHR == chr,'GEN'] <- sbPred + abs(min(sbPred))
        
      } else {
        
        pred[pred$CHR == chr,'GEN'] <- sbPred - min(sbPred)
      }
      
    } else if (what == 'gen2phys'){
      
      sbPred <- approx(x = fit$y, y = fit$x, xout = sbPred$POS)$y
      
      pred[pred$CHR == chr,'GEN'] <- round(sbPred)
      
    } else stop("The argument 'what' must be either 'phys2gen' or 'gen2phys'" )
    
  }
  
  return(pred)
}

# pred = cM_converter(map = 'MGC_marker_placements_for_genetic_map.map', pred = 'MGC_marker_placements_for_genetic_map.map', chrom_col = 1, pos_col = 4)
# write.table(pred, 'MGC_marker_placements_for_genetic_map_predicted.map', col.names=T, quote=F, append=F, row.names=F, sep='\t')

# pred = cM_converter(map = 'MGC_marker_placements_for_genetic_map.map', pred = 'MGC_imputed_maf05_markers.list')
# write.table(pred, 'MGC_imputed_maf05_markers.list', col.names=T, quote=F, append=F, row.names=F, sep='\t')


