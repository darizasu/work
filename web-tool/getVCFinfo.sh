#!/bin/bash

NGSEP=/data/software/NGSEPcore_3.0.2.jar
refIDs=/data/references/bean/v2.1/assembly/Pvulgaris_442_v2.0_seqnames.txt
VCFlocation=/bionas0/bean/populations/all_bean_summary/genotyping

samp1=$1
samp2=$2

echo -e 'chromosome\tposition\tref\talt\tsample1\tsample2\tdiffType\tcomparisons' > ${samp1}-${samp2}.tsv
java -Xmx10g -jar ${NGSEP} MergeVCF ${refIDs} ${VCFlocation}/${samp1}_bowtie2_NGSEP_gt.vcf.gz ${VCFlocation}/${samp2}_bowtie2_NGSEP_gt.vcf.gz 2>>${samp1}-${samp2}.log | java -jar ${NGSEP} FilterVCF -q 60 -minI 2 - 2>>${samp1}-${samp2}.log | vcf-query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' 2>>${samp1}-${samp2}.log | sed "s/Chr0//g" | sed "s/Chr//g" | grep -v 'scaffold' | python3 addZygocity.py >> ${samp1}-${samp2}.tsv 2>>${samp1}-${samp2}.log

################
test -s ${samp1}-${samp2}.log && echo 'Error: getVCFinfo failed at some point !!' && cat ${samp1}-${samp2}.log && rm ${samp1}-${samp2}.* && exit 1
################

Rscript --vanilla --restore CompareVCF_plot.R ${samp1}-${samp2}.tsv >> ${samp1}-${samp2}.log 2>&1

################
test -s ${samp1}-${samp2}.log && echo 'Error: CompareVCF_plot failed at some point !!' >> ${samp1}-${samp2}.log && cat ${samp1}-${samp2}.log && rm ${samp1}-${samp2}.tsv ${samp1}-${samp2}.log && exit 1
################

rm ${samp1}-${samp2}.tsv ${samp1}-${samp2}.log