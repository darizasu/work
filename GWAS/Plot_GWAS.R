#!/usr/bin/env Rscript


rm(list=ls())

args = commandArgs(trailingOnly=TRUE)

manhattan = function (x, chr = "CHR", bp = "BP", p = "P", snp = "SNP", col = c("gray10", "gray60"), chrlabs = NULL, suggestiveline = -log10(1e-05), genomewideline = -log10(5e-08), highlight = NULL, logp = TRUE, annotatePval = NULL, annotateTop = TRUE, ...){
  CHR = BP = P = index = NULL
  if (!(chr %in% names(x))) 
    stop(paste("Column", chr, "not found!"))
  if (!(bp %in% names(x))) 
    stop(paste("Column", bp, "not found!"))
  if (!(p %in% names(x))) 
    stop(paste("Column", p, "not found!"))
  if (!(snp %in% names(x))) 
    warning(paste("No SNP column found. OK unless you're trying to highlight."))
  if (!is.numeric(x[[chr]])) 
    stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
  if (!is.numeric(x[[bp]])) 
    stop(paste(bp, "column should be numeric."))
  if (!is.numeric(x[[p]])) 
    stop(paste(p, "column should be numeric."))
  d = data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]])
  if (!is.null(x[[snp]])) 
    d = transform(d, SNP = x[[snp]])
  d <- subset(d, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))
  d <- d[order(d$CHR, d$BP), ]
  if (logp) {
    d$logp <- -log10(d$P)
  }
  else {
    d$logp <- d$P
  }
  d$pos = NA
  d$index = NA
  ind = 0
  for (i in unique(d$CHR)) {
    ind = ind + 1
    d[d$CHR == i, ]$index = ind
  }
  nchr = length(unique(d$CHR))
  if (nchr == 1) {
    d$pos = d$BP
    ticks = floor(length(d$pos))/2 + 1
    xlabel = paste("Chromosome", unique(d$CHR), "position")
    labs = ticks
  }
  else {
    lastbase = 0
    ticks = NULL
    for (i in unique(d$index)) {
      if (i == 1) {
        d[d$index == i, ]$pos = d[d$index == i, ]$BP
      }
      else {
        lastbase = lastbase + tail(subset(d, index == 
                                            i - 1)$BP, 1)
        d[d$index == i, ]$pos = d[d$index == i, ]$BP + 
          lastbase
      }
      ticks = c(ticks, (min(d[d$index == i, ]$pos) + max(d[d$index == 
                                                             i, ]$pos))/2 + 1)
    }
    xlabel = "Chromosome"
    labs <- unique(d$CHR)
  }
  xmax = ceiling(max(d$pos) * 1.03)
  xmin = floor(max(d$pos) * -0.03)
  def_args <- list(xaxt = "n", bty = "n", xaxs = "i", yaxs = "i", 
                   las = 1, pch = 20, xlim = c(xmin, xmax), ylim = c(0, 
                                                                     ceiling(max(d$logp))), xlab = xlabel, ylab = expression(-log[10](italic(p))))
  dotargs <- list(...)
  do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% 
                                            names(dotargs)]))
  if (!is.null(chrlabs)) {
    if (is.character(chrlabs)) {
      if (length(chrlabs) == length(labs)) {
        labs <- chrlabs
      }
      else {
        warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
      }
    }
    else {
      warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
    }
  }
  if (nchr == 1) {
    axis(1, ...)
  }
  else {
    axis(1, at = ticks, labels = labs, ...)
  }
  col = rep(col, max(d$CHR))
  if (nchr == 1) {
    with(d, points(pos, logp, pch = 20, col = col[1], ...))
  }
  else {
    icol = 1
    for (i in unique(d$index)) {
      with(d[d$index == unique(d$index)[i], ], points(pos, 
                                                      logp, col = col[icol], pch = 20, ...))
      icol = icol + 1
    }
  }
  if (suggestiveline) 
    abline(h = suggestiveline, col = "blue")
  if (genomewideline) 
    abline(h = genomewideline, col = "red")
  if (!is.null(highlight)) {
    if (any(!(highlight %in% d$SNP))) 
      warning("You're trying to highlight SNPs that don't exist in your results.")
    d.highlight = d[which(d$SNP %in% highlight), ]
    with(d.highlight, points(pos, logp, col = "green3", 
                             pch = 20, ...))
  }
  if (!is.null(annotatePval)) {
    topHits = subset(d, P <= annotatePval)
    par(xpd = TRUE)
    if (annotateTop == FALSE) {
      with(subset(d, P <= annotatePval), textxy(pos, -log10(P), 
                                                offset = 0.625, labs = topHits$SNP, cex = 0.45), 
           ...)
    }
    else {
      topHits <- topHits[order(topHits$P), ]
      topSNPs <- NULL
      for (i in unique(topHits$CHR)) {
        chrSNPs <- topHits[topHits$CHR == i, ]
        topSNPs <- rbind(topSNPs, chrSNPs[1, ])
      }
      textxy(topSNPs$pos, -log10(topSNPs$P), offset = 0.625, 
             labs = topSNPs$SNP, cex = 0.5, ...)
    }
  }
  par(xpd = FALSE)
}

qq = function (pvector, ...){
  if (!is.numeric(pvector)) 
    stop("Input must be numeric.")
  pvector <- pvector[!is.na(pvector) & !is.nan(pvector) & 
                       !is.null(pvector) & is.finite(pvector) & pvector < 1 & 
                       pvector > 0]
  o = -log10(sort(pvector, decreasing = FALSE))
  e = -log10(ppoints(length(pvector)))
  def_args <- list(pch = 20, xlim = c(0, max(e)), ylim = c(0, 
                                                           max(o)), xlab = expression(Expected ~ ~-log[10](italic(p))), 
                   ylab = expression(Observed ~ ~-log[10](italic(p))))
  dotargs <- list(...)
  tryCatch(do.call("plot", c(list(x = e, y = o), def_args[!names(def_args) %in% 
                                                            names(dotargs)], dotargs)), warn = stop)
  abline(0, 1, col = "red")
}

plotManhattan = function(dataset, logfile, dataset_base, heritabilities = '/home/dariza/tassel_output/VEF_noMeso_mlm-PCA/heritabilities.txt'){
  
  setwd(getwd())
  
  ntaxa = system(paste("grep 'Number of Taxa' ",logfile," | grep -v 401", sep=''), intern = T)
  system(paste("sed -i 's/YDHA_2_05m.40g/YDHA_2_05m_40g/g'", dataset))
  system(paste("sed -i 's/YDHPL_2_10p.40g/YDHPL_2_10p_40g/g'", dataset))
  system(paste("sed '/SCAFFOLD/d'",dataset,">", paste(dataset,'.tmp',sep='')))
  myDataset = read.delim(paste(dataset,'.tmp',sep=''),colClasses = c(rep('factor',2),rep('integer',2),rep('NULL',2),'numeric',rep('NULL',11)))
  h2s = read.delim(heritabilities, row.names = 1)
  myDataset = myDataset[-1,]
  myDataset = myDataset[!is.na(myDataset[,3]),]
  row.names(myDataset) = seq(nrow(myDataset))
  colnames(myDataset) = c('TRAIT','SNP','CHR','BP','P')
  
  pdf(paste(dataset_base,'-plots.pdf',sep=''), width=12,height=4)
  
  for (trait in levels(myDataset[,1])){
    dataset_trait = myDataset[myDataset[,1] == trait,]
    
    if (all(is.na(dataset_trait$P))) next
    
    layout(matrix(c(1,1,2), nrow=1))
    manhattan(dataset_trait, ylim=c(0,10), main=paste(trait,dataset_base,sep=' - '), suggestiveline=F, genomewideline=F, cex = 0.5)
    bonf = 0.05/nrow(dataset_trait)
    abline(h=-log10(bonf), col='red')
    title(sub = ntaxa, adj = 0, line = 3, cex.sub = 0.85)
    name_2_h2s = gsub('_noMeso','',dataset_base)
    name_2_h2s = gsub('_byTassel','',name_2_h2s)
    title(sub = paste('h2 = ', h2s[name_2_h2s,trait],sep=''), adj = 1, line = 3, cex.sub = 0.85)
    qq(dataset_trait$P, main = paste(trait,dataset_base,sep=' - '))
    rm(dataset_trait)
    
  }
    garbage <- dev.off()
    system(paste("rm", paste(dataset,'.tmp',sep='')))
}

if (length(args)==0) {
  stop("At least three arguments must be supplied:\n
  \n  Usage: Rscript Plot_GWAS.R <TASSEL_output> <TASSEL_logFile> <title> \n
       Positional arguments:\n
       \tTASSEL_output     :  TASSEL output file with p-values for every marker. 
       \tTASSEL_logFile    :  TASSEL log file. 
       \ttitle             :  Base name used for the plots titles", call.=FALSE)
} else {
  plotManhattan(dataset = args[1], logfile = args[2], dataset_base = args[3])
}


