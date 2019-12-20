
f <- commandArgs(TRUE)[1]
i <- commandArgs(TRUE)[2]
x <- commandArgs(TRUE)[3]

s <- read.delim(text = paste0(head(readLines(paste0(f, '_bowtie2_readpos.stats')), -2), collapse='\n'), header = F)

f <- basename(f)

s$V6 = s$V2 / s$V4

o <- which((abs(s$V6 - median(s$V6)) / mad(s$V6)) > 3)

o <- R.utils::seqToIntervals(o)

if (nrow(o) == 1){
  
  cat(c(f, i, x, 0, o[nrow(o),1]), sep='\t', end='\n')
  
} else if (nrow(o) > 1){
  
  cat(c(f, i, x, o[1,2], o[nrow(o),1]), sep='\t', end='\n')
  
} else {
  
  cat(c(f, i, x, 0, 0), sep='\t', end='\n')
}

