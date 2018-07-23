rm(list=ls())

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) stop("At least one argument must be supplied (input file).n", call.=FALSE)

library(GENESIS)
library(GWASTools)
library(SNPRelate)

baseN = args[1]
traits = unlist(strsplit(args[2],','))
phen = args[3]

# Read genotypic data in plink BED format

if (! file.exists(paste(baseN,".gds",sep=''))){

  snpgdsBED2GDS(bed.fn    = paste(baseN,".bed",sep=''),
                bim.fn    = paste(baseN,".bim",sep=''),
                fam.fn    = paste(baseN,".fam",sep=''), 
                out.gdsfn = paste(baseN,".gds",sep=''))

}

geno = GdsGenotypeReader(filename = paste(baseN,".gds",sep=''))

genoData = GenotypeData(geno)


# Run PC-AiR and PC-Relate for PCA and Kinship matrix

if (! file.exists(paste(baseN,"_PC.RData",sep=''))){

  mypcair = pcair(genoData = genoData)
  save(mypcair, file=paste(baseN,"_PC.RData",sep=''))
}

load(paste(baseN,"_PC.RData",sep=''))

if (! file.exists(paste(baseN,"_kinship.RData",sep=''))){

  mypcrel = pcrelate(genoData = genoData, pcMat = mypcair$vectors[,1:5])
  save(mypcrel, file=paste(baseN,"_kinship.RData",sep=''))
}

load(paste(baseN,"_kinship.RData",sep=''))

mydat = data.frame(rows = mypcrel$sample.id,
                   scanID = mypcrel$sample.id,
                   pc1 = mypcair$vectors[,1], 
                   pc2 = mypcair$vectors[,2], 
                   pc3 = mypcair$vectors[,3], 
                   pc4 = mypcair$vectors[,4], 
                   pc5 = mypcair$vectors[,5],
                   row.names = 1)
phenoData = read.csv(phen, header = T, na.strings = 'NA', row.names=1)
cat('We got here')

mydat = cbind(mydat, phenoData[rownames(mydat),])

scanAnnot = ScanAnnotationDataFrame(mydat)

genoData = GenotypeData(geno, scanAnnot = scanAnnot)

myGRM = pcrelateMakeGRM(mypcrel)

for (trait in traits){
  
  baseN2 = paste(baseN,trait,sep='_')
  
  if (! file.exists(paste(baseN2,".results",sep=''))){
    
    nullmod = fitNullMM(scanData = scanAnnot,
                        outcome = trait,
                        covars = c('pc1','pc2','pc3','pc4','pc5'),
                        covMatList = myGRM, 
                        family = gaussian)
    
    assoc = assocTestMM(genoData = genoData, nullMMobj = nullmod, test = "Score")
    
    assoc$snpID = getVariable(geno, 'snp.rs.id')
    
    assoc$CHR = sapply(strsplit(x = assoc$snpID, split = ':'), "[[", 1)
    assoc$POS = sapply(strsplit(x = assoc$snpID, split = ':'), "[[", 2)
    assoc$REF = sapply(strsplit(x = assoc$snpID, split = ':'), "[[", 3)
    assoc$ALT = sapply(strsplit(x = assoc$snpID, split = ':'), "[[", 4)
    
    write.table(assoc, file=paste(baseN2,".results",sep=''), sep='\t', quote=F, row.names=F, col.names=T)
    
  }
  
  if (! file.exists(paste(baseN2,"_manhattan.pdf",sep=''))){
    
    initial = read.table(paste(baseN2,".results",sep=''), nrows = 100, comment.char='')
    classes = sapply(initial, class)
    assoc   = read.table(paste(baseN2,".results",sep=''), colClasses=classes, header=T, sep='\t')
    
    library(qqman)
    library(qvalue)
    
    pdf(paste(baseN2,"_manhattan.pdf",sep=''), width=12,height=4)
    
    if (all(is.na(assoc[,'Score.pval']))) next
    
    layout(matrix(c(1,1,2), nrow=1))
    
    ylim = c(0,ceiling(-log10(min(assoc[,'Score.pval'], na.rm=T))))
    
    manhattan(assoc, chr='chr', bp='POS', p='Score.pval', snp='snpID',
              suggestiveline=F, genomewideline=F, cex = 0.5, ylim=ylim, main=trait)
    
    bonf = 0.05/nrow(assoc)
    abline(h=-log10(bonf), col='red')
    
    fdr = qvalue(p = na.omit(assoc[,'Score.pval']))
    
    if (any(fdr$qvalues <= 0.05, na.rm=T)){
      
      fdr = max(fdr$pvalues[fdr$qvalues <= 0.05])
      abline(h=-log10(fdr), col='#26B14C')
      
    }
    
    qq(assoc[,'Score.pval'])
    
    garbage = dev.off()
    
  }
  
}
