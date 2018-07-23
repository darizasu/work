
rm(list=ls())

physGenMap = read.table('BARCBean_chip_marker_placements_for_genetic_map.map',
                        header=F,sep='\t', col.names = c('Chr','noidea','genetic','physical'))

physGenMap$Chr = factor(physGenMap$Chr)
fits = list()

for (chr in levels(physGenMap$Chr)){
  
  tmp_dataset = subset(physGenMap, Chr == chr)
  
  fits[[chr]] = lm(genetic ~ poly(physical, 10, raw = T), data = tmp_dataset)
  
  physGenMap[physGenMap$Chr == chr, 'predGenetic'] = sort(abs(predict(fits[[chr]])))
  physGenMap[physGenMap$Chr == chr, 'physical'] = sort(tmp_dataset$physical)
  
  cat(chr, cor(tmp_dataset$genetic, predict(fits[[chr]])), end='\n')
  # plot(tmp_dataset$genetic, predict(fits[[chr]]), main = chr)
  
}

physGenMap = physGenMap[with(physGenMap, order(Chr,predGenetic)),]

write.table(physGenMap[,c(1,2,5,4)],
            'BARCBean_chip_marker_placements_for_genetic_map_predicted.map',
            append=F, quote=F, sep='\t', row.names=F, col.names=F)
