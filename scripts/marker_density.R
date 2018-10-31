
marker_density <- function(markers, bin = 2.5e5, chromID = 'X.CHROM', posID = 'POS',
                           colorsRange = c('white','black'), chromLabs = NULL, 
                           pero_centro = NULL, plot_centro = F){
  
  # Object : Plot the density of markers from a genotypic matrix
  # Input  : markers . A table (may be gzip compressed) separated by tabs or commas.
  #                    This table contains the location of each marker in 
  #                    separate 'chromosome' and 'position' columns.
  # Input  : bin . Integer. Size of the color boxes in the plot in base-pairs.
  # Input  : chromID . Character. Name of the column in 'markers' with the chromosome info.
  # Input  : posID . Character. Name of the column in 'markers' with the position info.
  # Input  : colorsRange . Colors to interpolate; must be a valid argument to col2rgb()
  # Input  : chromLabs . Character vector with the labels to be used for the chromosomes.
  # Input  : pero_centro . A table separated by tabs or commas. It contains the coordinates 
  #                        for centromeric and pericentromeric regions in the genome.
  #                        It contains the columns 'pStart' and 'pEnd' for pericentromeres and/or
  #                        'cStart' and 'cEnd' for centromeres
  # Input  : plot_centro . Logical. Should the centromeric regions be drawn ?
  # Output : barplot object with the density of markers for each chromosome
  # Authors: darizasu
  #  Last update: October 31, 2018

  read_the_table <- function(ext, ...){
    if (ext == 'csv') return(read.csv(...)) else return(read.table(...))
  }

  if (tools::file_ext(markers) == 'gz') extCheck = tools::file_path_sans_ext(markers) else extCheck = markers

  if (tools::file_ext(extCheck) == 'csv') extCheck = 'csv' else extCheck = 'table'

  initial = read_the_table(ext = extCheck, file = markers, nrows = 10, comment.char = '##', header = T)
  classes = sapply(initial, class)
  index = grep( paste( c(chromID,posID), collapse = '|'), names(markers) )

  classes[-index] = "NULL"

  markers = read_the_table(ext = extCheck, file = markers, colClasses = classes, comment.char = '##', header = T)

  chroms <- unique(markers[,chromID])
  nchroms <- length(chroms)

  dens_matrix = matrix(nrow = 0, ncol = nchroms)
  cols_matrix = matrix(nrow = 0, ncol = nchroms)
  
  mxLen = c()

  for (i in chroms){

    Cmarkers = markers[ markers[,chromID] == i , posID ]
    mxLen = c(mxLen, max(Cmarkers))
    mx = ceiling( max(Cmarkers) / 1e6 ) * 1e6
    Cmarkers = cut( Cmarkers, b = seq(1, mx, bin), include.lowest = T)
    Cmarkers = table(Cmarkers)

    for (j in names(Cmarkers)){

      cols_matrix = rbind( cols_matrix, rep(NA, nchroms) )
      cols_matrix[ nrow(cols_matrix), i ] = Cmarkers[j]

      dens_matrix = rbind( dens_matrix, rep(0, nchroms) )
      dens_matrix[ nrow(dens_matrix), i ] = bin

    }
  }

  colfunc = colorRampPalette(colorsRange)
  ncols = na.omit( as.vector(cols_matrix) )
  cols = colfunc( max(ncols + 1))[ncols + 1]
  
  ylim = ceiling( max(markers[,posID], na.rm = T) / 1e6 ) * 1e6

  yTicks = seq(0, ylim, 1e7)
  yLabs  = paste(yTicks / 1e6, "Mb")
  
  if(is.null(chromLabs)) chromLabs = 1:nchroms
  

  par(mar=c(0.5, 5, 2, 1), fig = c(0,1,0,1))

  bp = barplot(dens_matrix, horiz = F, las = 2, xaxt = 'n', yaxt = 'n',
               border = NA, space = 2, cex.names = 1, ylab = 'Position',
               col = cols, ylim = c(ylim, 0))

  axis(3, labels = chromLabs, at = bp, las=1, cex.axis=0.75, tick=F, line=-.5, font=2)

  axis(2, at = yTicks, labels = yLabs, las=2, cex.axis=0.75)
  
  xSegm = c(rbind(seq(3,33,3) - 1, seq(3,33,3)))
  ySegm = c(rbind(mxLen,mxLen))
  segments(x0 = xSegm, y0 = 0, x1 = xSegm, y1 = ySegm)
  
  xSegm = c(rbind(seq(3,33,3),seq(3,33,3)))
  ySegm = c(rbind(rep(0, nchroms), mxLen))
  segments(x0 = xSegm - 1, y0 = ySegm, x1 = xSegm, y1 = ySegm)
  
  if(!is.null(pero_centro)){

    if (tools::file_ext(pero_centro) == 'gz') extCheck = tools::file_path_sans_ext(pero_centro) else extCheck = pero_centro
    if (tools::file_ext(extCheck) == 'csv') extCheck = 'csv' else extCheck = 'table'

    segTable = read_the_table(ext = extCheck, file = pero_centro, comment.char = '##', header = T)
    
    pericent = c('pStart','pEnd') %in% colnames(segTable)
    
    if(sum(pericent)){
      
      ySegm = c(rbind(segTable$pStart, segTable$pEnd))
      segments(x0 = xSegm - 1, y0 = ySegm, x1 = xSegm, y1 = ySegm)
      
    } else warning("The pero_centro table does not contain the columns 'pStart' or 'pEnd'")
    
    cent = c('cStart','cEnd') %in% colnames(segTable)
    
    if(sum(cent) * plot_centro){
      
      ySegm = c(rbind(segTable$cStart, segTable$cEnd))
      
      segments(x0 = xSegm - 1, y0 = ySegm, x1 = xSegm,
               y1 = c(rbind(ySegm[c(F,T)],ySegm[c(T,F)])))
      
      xSegm = c(rbind(seq(3,33,3) - 1, seq(3,33,3)))
      cStart = c(rbind(segTable$cStart,segTable$cStart)) + 3e5
      cEnd = c(rbind(segTable$cEnd,segTable$cEnd)) - 3e5
      segments(x0 = xSegm, y0 = cStart, x1 = xSegm, y1 = cEnd, col = 'white')
      
    } else warning("The centromeres will not be drawn")
  }

  par(fig = c(.4, .6, .1, .2), new = T)

  image(x = 1:length(ncols),z = t(t(1:length(ncols))),
        ylim = c(0,1), axes = F, xlab = "",
        col = colfunc( max(ncols + 1)))

  axis(1, at = c(1,length(ncols)), labels = c(0,length(ncols)), tick = T, cex.axis=0.6)
  
}

# png('introgressions_all.png', width = 8, height = 6, units = 'in', res = 300)
# 
# marker_density(markers = 'MGC_marker_positions_GBS.txt', bin = 2.5e5, chromID = 'X.CHROM', posID = 'POS',
#                chromLabs = c('Pv01','Pv02','Pv03','Pv04','Pv05','Pv06','Pv07','Pv08','Pv09','Pv10','Pv11'),
#                pero_centro = 'Pvulgaris_centromeres.txt')
# dev.off()
