rm(list = ls())

pks <- c('dplyr','argparse')

pks <- suppressPackageStartupMessages(sapply(pks, require, character.only = TRUE))

if (any(!pks)) stop('The package(s) ', paste(names(pks)[!pks]), ' is/are not available for load.')

pr <- ArgumentParser(description = "This script calculates heterozygosity rates per sample and per marker from a VCF file. It must be used after parsing the VCF like: bcftools query -H -f '%ID[\\t%GT]\\n' myVCF.vcf.gz",
                     formatter_class= 'argparse.RawTextHelpFormatter')

pr$add_argument("-i", type = "character", metavar = "file", default = "stdin",
                help = "Input GT matrix in the format described above")

pr$add_argument("-o", type = "character", metavar = "file",
                help = "Output PDF file")

pr$add_argument("--save-table", action="store_true", default=FALSE,
                help="Save heterozygosity rates in CSV files [Default: Do not save tables]")

ar = pr$parse_args()

if ( any(sapply(ar, is.null)) ){
  
  # Display the help message and exit
  system( paste( paste(commandArgs()[1:4], collapse = " "), 
                 '--args --help', collapse = "" ) )
  stop('One or more arguments are not valid. Check usage for more details.')
}

input <- ar$i
output <- ar$o
sv <- ar$save_table

H <- 
  read.delim(file = input, header = T, na.strings = "./.", check.names = F)

sam <-
  colnames(H) %>% 
  gsub("\\[.*?\\]","", .) %>% 
  gsub(":GT", "", .)

sam <- sam[-1]
colnames(H) <- paste0('V', 1:ncol(H))

if (n_distinct(H$V1) > 5){
  
  ids <- H[,1]
  H   <- H[,-1]
}

H <- 
  mutate_all(H, (list( ~ case_when( . == '0/0' ~ -1,
                                     . == '0|0' ~ -1,
                                     . == '0/1' ~  0,
                                     . == '0|1' ~  0,
                                     . == '1|0' ~  0,
                                     . == '1/1' ~  1,
                                     . == '1|1' ~  1) ))) %>% 
  as.matrix()

H <- 1 - abs(H)

H.snp <- apply(H, 1, mean, na.rm = TRUE)
H.sam <- apply(H, 2, mean, na.rm = TRUE)

pdf("Heterozygosity.pdf", width = 10, height = 6)

par(mfrow=c(1,2), mar=c(5,5,1,1)+0.1)

hist(H.sam,col = "gray", main = "", xlab = "Heterozygosity of individuals",
     ylab = paste0("Frequency (out of ", length(H.sam)," individuals)"))

hist(H.snp,col = "gray", main = "", xlab = "Heterozygosity of markers",
     ylab = paste0("Frequency (out of ", length(H.snp)," markers)"))

dev.off()

if (sv){
  
  write.csv(x = data.frame(Sample = sam, Heterozygosity = H.sam),
            file = 'Het.samples.csv', row.names = F, quote = F)
  
  if (exists("ids")){
    
    write.csv(x = data.frame(SNP = ids, Heterozygosity = H.snp),
              file = 'Het.SNPs.csv', row.names = F, quote = F)
  }
}

