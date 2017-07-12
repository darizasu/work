
# Stacked barplot tile for gv2.1 ------------------------------------------

rm(list = ls())

setwd("F:/Bean_Introgressions_paper")
# introgressions = read.table('introgressions_v2.1.txt', header = T, sep = '\t')
introgressions = read.table('introgressions_v2.1_1Mb.txt', header = T, sep = '\t')

introgressions$sample = factor(introgressions$sample, levels = c("AND_696","G19833","G24705","G4627","G5686","Gasilida","Ibaddo","ICA_Quimbaya","Kablanketi","Midas","Mshindi","NABE_12C","NABE_13","NABE_14","ALB_213","G10474","G2333","G750","Hawassa_Dume","INB_827","INB_841","Mexico_54","MIB_778","Red_Wolayta","SCR_2","SCR_9","SEA_5","SEN_56","SXB_412","TARS-VR-7S","VAX_1","VAX_3","VAX_4","VAX_6"))

chrom_lengths = data.frame(Pv01=51433939,Pv02=49670989,Pv03=53438756,Pv04=48048378,Pv05=40923498,Pv06=31236378,Pv07=40041001,Pv08=63048260,Pv09=38250102,Pv10=44302882,Pv11=53580169)

samp_background = c(rep('Andean', 14), rep('Mesoamerican', 20))
names(samp_background) = levels(introgressions$sample)


myList = list()

# pdf("Introgressions_BeanWGS_v2.1.pdf",15,7)

for(chrom in levels(introgressions$chromosome)){
# for(chrom in c('Pv08')){
  
  # png(paste(chrom, '.png', sep=''), width=30,height=9, units='in', res=300)
  
  myList[[chrom]] = matrix(ncol = nlevels(introgressions$sample))
  col = paste(chrom,'_col', sep = '')
  # myList[[col]] = vector('character')
  myList[[col]] = matrix(ncol = nlevels(introgressions$sample))
  colnames(myList[[chrom]]) = levels(introgressions$sample)
  colnames(myList[[col]]) = levels(introgressions$sample)
  
  for(sample in levels(introgressions$sample)){
    
    tmp_subset = introgressions[introgressions[,'chromosome'] == chrom & introgressions[,'sample'] == sample,]

    if(nrow(tmp_subset) == 0){
      
      back_group = as.character(samp_background[sample])
      myList[[chrom]] = rbind(myList[[chrom]], rep(0,nlevels(introgressions$sample)))
      myList[[chrom]][nrow(myList[[chrom]]),sample] = chrom_lengths[,chrom]
      # myList[[col]] = c(myList[[col]], back_group)
      myList[[col]] = rbind(myList[[col]], rep(0,nlevels(introgressions$sample)))
      myList[[col]][nrow(myList[[col]]),sample] = back_group
      
    } else {
      
      tmp_subset = tmp_subset[order(tmp_subset$start),]
      tmp_subset$back_group = as.character(tmp_subset$back_group)
      tmp_subset$hap_group = as.character(tmp_subset$hap_group)
      back_group = as.character(unique(tmp_subset$back_group))
      positions = sort(c(tmp_subset$start, tmp_subset$end, 0, chrom_lengths[,chrom]))
      # myList[[col]] = c(myList[[col]],rbind(tmp_subset$back_group,tmp_subset$hap_group),back_group)
      group = c(rbind(tmp_subset$back_group,tmp_subset$hap_group),back_group)
      if (any(positions > chrom_lengths[,chrom])){
        positions = positions[positions<=chrom_lengths[,chrom]]
        group = c(rbind(tmp_subset$back_group,tmp_subset$hap_group))
      }

      for(i in 2:length(positions)){
        
        myList[[chrom]] = rbind(myList[[chrom]], rep(0,nlevels(introgressions$sample)))
        myList[[chrom]][nrow(myList[[chrom]]),sample] = positions[i]-positions[i-1]
        myList[[col]] = rbind(myList[[col]], rep(0,nlevels(introgressions$sample)))
        myList[[col]][nrow(myList[[col]]),sample] = group[i-1]
        
      }
    }
  }
  
  myList[[chrom]] = myList[[chrom]][2:nrow(myList[[chrom]]),]
  myList[[col]] = myList[[col]][2:nrow(myList[[col]]),]
  myList[[col]] = as.vector(myList[[col]])
  myList[[col]] = myList[[col]][myList[[col]] != 0]
  myList[[col]][myList[[col]] == 'Andean'] = '#4285F4'
  myList[[col]][myList[[col]] == 'Mesoamerican'] = '#FFCCC9'
  myList[[col]][myList[[col]] == 'Acutifolius'] = '#34A853'
  myList[[col]][myList[[col]] == 'Coccineus'] = '#E1A904'
  
  # par(mar=c(0.5, 4, 6.5, 2))
  par(mar=c(0.5, 4, 0.5, 1.5))
  
  bp = barplot(myList[[chrom]], col=myList[[col]], horiz=F, las=2, border=F, space=0.1, ylim = c(chrom_lengths[,chrom],0), xaxt='n', yaxt='n',cex.names=1)
  # mtext(chrom, side = 4, font = 2, cex = 1.5)
  
  corners = par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
  par(xpd = TRUE) #Draw outside plot area
  text(x = corners[2]+.5, y = mean(corners[3:4]), chrom, srt = 270, font = 2, cex = 1.5)
  
  axis(3, labels = F, at = bp, las=2, cex.axis=1, tick=F, line=-0.8, font=2)
  # axis(3, labels = gsub('_', ' ', levels(introgressions$sample)), at = bp, las=2, cex.axis=1, tick=F, line=-0.8, font=2)
  
  index = chrom_lengths[,chrom] %/% 1e7
  myLabels = c('0 Mb', '10 Mb','20 Mb','30 Mb','40 Mb','50 Mb','60 Mb')
  
  axis(2, at = seq(0, index*1e7, 1e7), labels = myLabels[0:index+1], las=2)
  # axis(2, at = seq(0, index*1e7, 1e7), labels = F, las=2)
# dev.off()
}

# dev.off()

# Stacked barplot tile for gv2.1 with adjusted bar length -----------------

rm(list = ls())

setwd("F:/Bean_Introgressions_paper")
introgressions = read.table('introgressions_v2.1_1Mb.txt', header = T, sep = '\t')

introgressions$sample = factor(introgressions$sample, levels = c("AND_696","G19833","G24705","G4627","G5686","Gasilida","Ibaddo","ICA_Quimbaya","Kablanketi","Midas","Mshindi","NABE_12C","NABE_13","NABE_14","ALB_213","G10474","G2333","G750","Hawassa_Dume","INB_827","INB_841","Mexico_54","MIB_778","Red_Wolayta","SCR_2","SCR_9","SEA_5","SEN_56","SXB_412","TARS-VR-7S","VAX_1","VAX_3","VAX_4","VAX_6"))

chrom_lengths = data.frame(Pv01=51433939,Pv02=49670989,Pv03=53438756,Pv04=48048378,Pv05=40923498,Pv06=31236378,Pv07=40041001,Pv08=63048260,Pv09=38250102,Pv10=44302882,Pv11=53580169)

samp_background = c(rep('Andean', 14), rep('Mesoamerican', 20))
names(samp_background) = levels(introgressions$sample)

whitespace = 1.5e6

plot_per_Chr = function(table, chromLengths, sampBackground, chroms2plot){
  
  table = table[table$chromosome %in% chroms2plot,]
  for (column in colnames(table)){
    if (is.factor(table[,column])) table[,column] = factor(table[,column])
  }
  
  
  coord_matrix = matrix(ncol = nlevels(table$sample))
  color_matrix = matrix(ncol = nlevels(table$sample))
  colnames(coord_matrix) = levels(table$sample)
  colnames(color_matrix) = levels(table$sample)
  
  for(sample in levels(table$sample)){
    
    for(chrom in levels(table$chromosome)){
      
      tmp_subset = table[table[,'chromosome'] == chrom & table[,'sample'] == sample,]
      
      if(nrow(tmp_subset) == 0){
        
        back_group = as.character(sampBackground[sample])
        coord_matrix = rbind(coord_matrix, rep(0,nlevels(table$sample)))
        coord_matrix = rbind(coord_matrix, rep(0,nlevels(table$sample)))
        coord_matrix[nrow(coord_matrix)-1,sample] = chromLengths[,chrom]
        coord_matrix[nrow(coord_matrix),sample] = whitespace
        
        color_matrix = rbind(color_matrix, rep(0,nlevels(table$sample)))
        color_matrix = rbind(color_matrix, rep(0,nlevels(table$sample)))
        color_matrix[nrow(color_matrix)-1,sample] = back_group
        color_matrix[nrow(color_matrix),sample] = 'transparent'
        
      } else {
        
        tmp_subset = tmp_subset[order(tmp_subset$start),]
        tmp_subset$back_group = as.character(tmp_subset$back_group)
        tmp_subset$hap_group = as.character(tmp_subset$hap_group)
        back_group = as.character(unique(tmp_subset$back_group))
        positions = sort(c(tmp_subset$start, tmp_subset$end, 0, chromLengths[,chrom]))
        group = c(rbind(tmp_subset$back_group,tmp_subset$hap_group),back_group)
        if (any(positions > chromLengths[,chrom])){
          positions = positions[positions<=chromLengths[,chrom]]
          group = c(rbind(tmp_subset$back_group,tmp_subset$hap_group))
        }
        
        for(i in 2:length(positions)){
          
          coord_matrix = rbind(coord_matrix, rep(0,nlevels(table$sample)))
          coord_matrix[nrow(coord_matrix),sample] = positions[i]-positions[i-1]
          
          color_matrix = rbind(color_matrix, rep(0,nlevels(table$sample)))
          color_matrix[nrow(color_matrix),sample] = group[i-1]
          
        }
        
        coord_matrix = rbind(coord_matrix, rep(0,nlevels(table$sample)))
        coord_matrix[nrow(coord_matrix),sample] = whitespace
        
        color_matrix = rbind(color_matrix, rep(0,nlevels(table$sample)))
        color_matrix[nrow(color_matrix),sample] = 'transparent'
        
      }
    }
    
  }
  coord_matrix = coord_matrix[2:nrow(coord_matrix),]
  color_matrix = color_matrix[2:nrow(color_matrix),]
  color_matrix = as.vector(color_matrix)
  color_matrix = color_matrix[color_matrix != 0]
  color_matrix[color_matrix == 'Andean'] = '#4285F4'
  color_matrix[color_matrix == 'Mesoamerican'] = '#FFCCC9'
  color_matrix[color_matrix == 'Acutifolius'] = '#34A853'
  # color_matrix[color_matrix == 'Coccineus'] = '#E1A904'
  color_matrix[color_matrix == 'Coccineus'] = '#EA4335'
  
  par(mar=c(0.5, 5, 5.5, 1))
  # par(mar=c(0.5, 4, 6.5, 2))
  # par(mar=c(0.5, 4, 0.5, 1.5))
  
  bp = barplot(coord_matrix, col=color_matrix, horiz=F, las=2, border=F, space=0.1, ylim = c(sum(chromLengths, whitespace*10),0), xaxt='n', yaxt='n',cex.names=1)
  # mtext(chrom, side = 4, font = 2, cex = 1.5)
  
  axis(3, labels = gsub('_', ' ', levels(table$sample)), at = bp, las=2, cex.axis=0.65, tick=F, line=-0.8, font=2)
  
  index = c(5,5,5,5,4,3,4,6,4,4,6)
  myLabels = c('0 Mb', '10 Mb','20 Mb','30 Mb','40 Mb','50 Mb','60 Mb')
  allLabels = vector('character')
  for (i in index) {allLabels = c(allLabels, myLabels[1:i])}
  start_point = 0
  for (i in 2:ncol(chromLengths)) {start_point = c(start_point, start_point[i-1] + whitespace + chromLengths[,i-1])}
  allLabels_pos = vector('numeric')
  for (i in 1:length(index)) {
    for (j in 1:index[i]){
      allLabels_pos = c(allLabels_pos, ((j-1)*1e7) + start_point[i])
    }
  }
  
  axis(2, at = allLabels_pos, labels = allLabels, las=2, cex.axis=0.75)
  # axis(2, at = seq(0, index*1e7, 1e7), labels = F, las=2)
  
  corners = par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
  par(xpd = TRUE) #Draw outside plot area
  # text(x = corners[2]+.5, y = mean(corners[3:4]), chrom, srt = 270, font = 2, cex = 1.5)
  text(x = -8.6, y = as.vector((chrom_lengths / 2) + start_point), levels(table$chromosome), srt = 90, font = 2, cex = 1.1)
  bp
  
}


png(paste('introgressions_all', '.png', sep=''), width=6,height=16, units='in', res=1000)
# layout(matrix(1:2, ncol=2))
plot_per_Chr(table = introgressions, chromLengths = chrom_lengths, sampBackground = samp_background, names(chrom_lengths))
# plot_per_Chr(table = introgressions, chromLengths = chrom_lengths, sampBackground = samp_background, names(chrom_lengths))

dev.off()

rm(list = ls())

# SNP density -------------------------------------------------------------

rm(list = ls())

setwd("F:/Bean_Introgressions_paper")

sampleTable <- read.delim('SNP_density_positions.txt')#, colClasses = myClasses)

start_centr = c(6.8e6,4.5e6,6e6,8e6,4e6,0e6,9.8e6,9.8e6,1.5e6,5.2e6,9.8e6)
end_centr = c(38e6,25.5e6,29.5e6,39.5e6,33.8e6,15e6,37.5e6,48e6,5.8e6,34e6,43e6)

png("SNP_density.png",width = 15,height = 7, units = 'in', res = 600)
LSD::heatscatter(sampleTable[,1], sampleTable[,2], colpal = 'standardheat', cexplot = 3, pch = '_', xlab = 'Chromosome',
                 ylab = 'Position', main = "SNP_density", rev = T, ylim = c(6.5e7,0))

chroms = c(rbind(seq(1,11,1),seq(1,11,1)))
from.x = chroms - 0.2
to.x = chroms + 0.2
y = c(rbind(c(6.8e6,4.5e6,6e6,8e6,4e6,0e6,9.8e6,9.8e6,1.5e6,5.2e6,9.8e6),
            c(38e6,25.5e6,29.5e6,39.5e6,33.8e6,15e6,37.5e6,48e6,5.8e6,34e6,43e6)))

segments(x0 = from.x, y0 = y, x1 = to.x, y1 = y)

dev.off()

# png("barcolor_scale.png",width = 15,height = 7, units = 'in', res = 600)
# LSD::disco()
# dev.off()
