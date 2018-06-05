#!/home/dariza/bin/Rscript

rm(list=ls())

lop = c("argparse","tools")
new.pckgs = lop[!(lop %in% rownames(installed.packages()))]

if (length(new.pckgs)){
  cat("\nWarning message:\nThis script will try to install the packages",
  	  new.pckgs,"and all its dependencies. Otherwise abort this execution (Crtl + C)\n")
  Sys.sleep(10)
  install.packages(pkgs=new.pckgs, repos="http://cran.r-project.org", dependencies = T)
}

library(argparse)

parser = ArgumentParser()

parser$add_argument("-f", type="character", metavar='file',
                    help="File with GWAS results. It must have a header line. ")
parser$add_argument("-p", type='character', metavar='pvalues',
                    help="Name of the column that contains the p-values from the GWAS. This column must contain numerical values. Missing values should be specified as NA")

args = parser$parse_args()


if (any(sapply(args, is.null))){
  
  stop('One or more arguments are not valid. Check usage for more details.')
  
}

read_the_table = function(ext, ...){
  if (ext == 'csv') return(read.csv(...)) else return(read.table(...))
}

myFile = args$f
myP = args$p

if (tools::file_ext(myFile) == 'gz') extCheck = tools::file_path_sans_ext(myFile) else extCheck = myFile

if (tools::file_ext(extCheck) == 'csv') extCheck = 'csv' else extCheck = 'table'

initial = read_the_table(ext=extCheck, file=myFile, nrows=100, comment.char='', header=T)
classes = sapply(initial, class)
index   = grep(myP, names(classes))

classes[-index] = "NULL"
myP     = read_the_table(ext=extCheck, file=myFile, colClasses=classes, header=T)[,1]
myP     = na.omit(myP)
myP     = sort(myP, decreasing=F)

myE     = ppoints(length(myP))

MSD  = mean( sqrt( (myP - myE) ^ 2 ) )
RMSE = sqrt(  sum( (myP - myE) ^ 2 ) / length(myP) )

cat('\nMSD is' , MSD, '\n')
cat('RMSE is', RMSE, '\n\n')
