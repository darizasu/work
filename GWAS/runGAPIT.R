library(multtest)
library(gplots)
library(LDheatmap)
library(genetics)
library(ape)
library(EMMREML)
library(compiler) #this library is already installed in R
library("scatterplot3d")
# library(data.table)

source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("http://zzlab.net/GAPIT/emma.txt")

setwd("/bioinfo1/projects/bean/MGC/GWAS/GAPIT/MGC_test_2018")

myY = read.csv('/bioinfo1/projects/bean/MGC/GWAS/GAPIT/BLUPs_MGC_all.csv', header = T)
myY = myY[,-2]

initial = read.table('/bioinfo1/projects/bean/MGC/GWAS/GAPIT/MGC_imputed_fi_fir_maf06_oh06_hmp.txt.gz', nrows = 100, comment.char='')
classes = sapply(initial, class)
myG     = read.table('/bioinfo1/projects/bean/MGC/GWAS/GAPIT/MGC_imputed_fi_fir_maf06_oh06_hmp.txt.gz', colClasses=classes, header=F, comment.char='')

myCV = read.csv('GAPIT.PCA.csv', header = T)
myKI = read.csv('GAPIT.Kin.VanRaden.csv', header = F)

mygapit = GAPIT(
  Y=myY,
  G=myG,
  CV=myCV,
  KI=myKI
  # PCA.total=5,
  # kinship.algorithm = "EMMA",
  # file.fragment = 128,
  # group.from=426,
  # group.to=426,
  # group.by=1
)

