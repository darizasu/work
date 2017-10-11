args = commandArgs(trailingOnly=TRUE)
vargs <- strsplit(args, ",")

phen   = vargs[[1]][1]
phen2  = vargs[[1]][2]

str(phen)
!is.na(phen2)

length(vargs[[1]])

if (!is.na(phen2) && length(vargs[[1]]) == 2) print ('access')
