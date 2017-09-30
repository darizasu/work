#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

plotDiffs <- function(datasetFullpath, title) {

  myClasses <- c("numeric","numeric","character","character","factor")
  sampleTable <- read.delim(datasetFullpath, colClasses = myClasses)

  sampleTable[sampleTable[,5] == 'HeteroDiff', 1] <- sampleTable[sampleTable[,5] == 'HeteroDiff', 1] - 0.15
  sampleTable[sampleTable[,5] == 'HomoDiff', 1] <- sampleTable[sampleTable[,5] == 'HomoDiff', 1] + 0.15
  
  noHetero <- length(sampleTable[sampleTable[,5] == 'HeteroDiff', 1])
  noHomo <- length(sampleTable[sampleTable[,5] == 'HomoDiff', 1])
  noTotal <- length(sampleTable[,1])
  Hetero <- round((noHetero * 100) / noTotal, 2)
  Homo <- round((noHomo * 100) / noTotal, 2)

  pdf(paste(title, ".pdf", sep=''),15,7)
  LSD::heatscatter(sampleTable[,1], sampleTable[,2], colpal = 'bl2rd', cexplot = 0.8, pch = '_', xlab = 'Chromosome',
              ylab = 'Position', main = title, panel.first = abline(h = seq(5e6, 6e7, 5e6), col = '#E7E7E7'))
  title(sub = paste(Homo, '% Homozygous Differences', sep = ''), adj = 1, line = 2, cex.sub = 0.85)
  title(sub = paste(Hetero, '% Heterozygous Differences', sep = ''), adj = 1, line = 3, cex.sub = 0.85)
  title(sub = paste(noTotal, ' genotype calls in common', sep = ''), adj = 0, line = 2, cex.sub = 0.85)
  garbage <- dev.off()
}

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else {
  plotDiffs(datasetFullpath = args[1], title = args[2])
}
