#!/bin/bash

VCFlocation=/bionas1/bean/populations/all_bean_summary/genotyping

samp1=$1
samp2=$2
Q=60

################
test ! -s ${VCFlocation}/${samp1}_bowtie2_NGSEP_gt.vcf.gz && echo 'Error: The file '${VCFlocation}/${samp1}_bowtie2_NGSEP_gt.vcf.gz' does not exist' && exit 1
test ! -s ${VCFlocation}/${samp2}_bowtie2_NGSEP_gt.vcf.gz && echo 'Error: The file '${VCFlocation}/${samp2}_bowtie2_NGSEP_gt.vcf.gz' does not exist' && exit 1
################

python3 Compare2VCFs.py -q ${Q} ${VCFlocation}/${samp1}_bowtie2_NGSEP_gt.vcf.gz ${VCFlocation}/${samp2}_bowtie2_NGSEP_gt.vcf.gz > ${samp1}-${samp2}.tsv 2>${samp1}-${samp2}.log

################
test -s ${samp1}-${samp2}.log && echo 'Error: getVCFinfo failed at some point !!' && cat ${samp1}-${samp2}.log && rm ${samp1}-${samp2}.* && exit 1
################

Rscript --vanilla --restore CompareVCF_plot.R ${samp1}-${samp2}.tsv ${samp1}-${samp2} >> ${samp1}-${samp2}.log 2>&1

################
test -s ${samp1}-${samp2}.log && echo 'Error: CompareVCF_plot failed at some point !!' >> ${samp1}-${samp2}.log && cat ${samp1}-${samp2}.log && rm ${samp1}-${samp2}.tsv ${samp1}-${samp2}.log && exit 1
################

rm ${samp1}-${samp2}.tsv ${samp1}-${samp2}.log