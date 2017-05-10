#!/bin/bash

  # Key variables to specify

WD=/bioinfo2/projects/bean_tmp/organizing_data_tmp/populations/AXM
popName=AXM
quality=(60); # Specify the quality thresholds for FilterVCF
numThreads=20

# The following file should be located at ${WD}/genotyping, otherwise its path must be specified.
# This file contains four tab-separated columns:
# /path/to/sample[tab]sample_name[tab]ignore5[tab]ignore3
# '/path/to/sample' must contain a VCF and a BAM file for each sample on the same specified directory
# 'sample_name' is the name the specified sample will take in the final VCF file.
# You must name every sample as (e.g.): samplename-p01F12, the '-p' is MANDATORY after 'samplename', and 01F12 may indicate the plate and well that originated that sample.
# ignore5 is to ignore this many base pairs from the 5' end in NGSEP FindVariants
# ignore3 is to ignore this many base pairs from the 3' end in NGSEP FindVariants. See the example below:
# /bioinfo2/projects/GBSplates/01/mapping/ALB_213	ALB_213-p01H10	4	10
# This file may contain comment lines with '#'
samples2compare=/bioinfo2/projects/bean_tmp/organizing_data_tmp/populations/AXM/samples2merge_AXM.txt; # These are all the samples from your population

referenceIDs=/data/references/bean/v2.1/assembly/Pvulgaris_442_v2.0_seqnames.txt
knownVariants=/data/references/bean/v2.1/variants/GbS_p1-20.vcf

  # Path to Software used

NGSEP=/data/software/NGSEPcore_3.0.2.jar
samtoolsExec=/data/software/samtools/samtools-1.2

  # Reference genome files

REF=/data/references/bean/v2.1/bowtie2/Pvulgaris_442_v2.0.fa
STRs=/data/references/bean/v2.1/strs/Pvulgaris_v2_strs.list


#### ------------------------------------------------------------------------------------------------------------------------------- ####
#### ------------------------------------------------------------------------------------------------------------------------------- ####


# Create a single file containing the samples to be included in every VCF.

cd ${WD}/genotyping

samples=(`grep -v '^#' ${samples2compare} | cut -f 2 | sed 's/-p.*//' | tr '\n' ' '`)
sampleIDs=(`grep -v '^#' ${samples2compare} | cut -f 2 | sed 's/-p.*//' |  sort -d | uniq -d | tr '\n' ' '`)
locations=(`grep -v '^#' ${samples2compare} | cut -f 1 | tr '\n' ' '`)
list=(`grep -v '^#' ${samples2compare} | cut -f 2 | tr '\n' ' '`)

echo 'This run contains the following samples:'; echo ${sampleIDs[@]}; echo ''

for index in ${!samples[@]}
do

  if [[ " ${sampleIDs[@]} " =~ " ${samples[${index}]} " ]]
  then

    index2=`expr ${index} + 1`
    grep -v '^#' ${samples2compare} | sed -n -e ${index2}p >> ${samples[${index}]}.tmp

  else

    :

  fi

done

myTMPs=(`ls *tmp`)

for index in ${!myTMPs[@]}; do samplesPerList=`expr ${index} % ${numThreads}`; echo ${myTMPs[${index}]} >> TMPsList_${samplesPerList}.tmp;  done


#### ------------------------------------------------------------------------------------------------------------------------------- ####
#### ------------------------------------------------------------------------------------------------------------------------------- ####


cd ${WD}/genotyping


#### -------- Find Variants ------- ####


for file in ${WD}/genotyping/TMPsList*.tmp
do

  (

  TMPs=(`cat ${file} | tr '\n' ' '`)

  for sample in ${TMPs[@]}
  do

    sleep 5
    locations=(`grep -v '^#' ${sample} | cut -f 1 | tr '\n' ' '`)
    list=(`grep -v '^#' ${sample} | cut -f 2 | tr '\n' ' '`)
    Ifive=(`grep -v '^#' ${sample} | cut -f 3 | tr '\n' ' '`)
    i5=`echo ${Ifive[@]} | tr ' ' '\n' | sort -nr | head -1`
    Ithree=(`grep -v '^#' ${sample} | cut -f 4 | tr '\n' ' '`)
    i3=`echo ${Ithree[@]} | tr ' ' '\n' | sort -nr | head -1`
    sampName=`grep -v '^#' ${sample} | cut -f 2 | sed 's/-p.*//' | uniq`

    echo 'Starting find variants on '${sampName}' files' $(date)

    mkdir ${sampName}

    for index in ${!list[@]}
    do

      ln -s ${locations[${index}]}_bowtie2_sorted.bam ./${sampName}/${list[${index}]}_bowtie2_sorted.bam

      echo $(date) 'Genotyping population variants in '${list[${index}]}

      java -Xmx3g -jar ${NGSEP} FindVariants -h 0.0001 -maxBaseQS 30 -minQuality 0 -noRep -noRD -noRP -maxAlnsPerStartPos 100 -ignore5 ${Ifive[${index}]} -ignore3 ${Ithree[${index}]} -sampleId ${list[${index}]} -knownVariants ${knownVariants} ${REF} ./${sampName}/${list[${index}]}_bowtie2_sorted.bam ./${sampName}/${list[${index}]}_bowtie2_NGSEP_gt >& ./${sampName}/${list[${index}]}_bowtie2_NGSEP_gt.log

      ################
      if [[ ! `tail -1 ./${sampName}/${list[${index}]}_bowtie2_NGSEP_gt.log` == *Completed* ]]; then echo "Error: Find variants in "${list[${index}]}" failed at some point !!"; exit 1; fi 
      ################

    done

    echo 'Genotyping population variants on '${sampName}' files seems to be completed successfully' $(date)


  #### --------- Merge VCF ---------- ####

    echo 'Starting MergeVCF '${sampName}' files' $(date)

    java -Xmx10g -jar ${NGSEP} MergeVCF ${referenceIDs} ${WD}/genotyping/${sampName}/${sampName}*_gt.vcf > ${WD}/genotyping/${sampName}/${sampName}.vcf

    ################
    test ! -s ${WD}/genotyping/${sampName}/${sampName}.vcf && echo "Error: MergeVCF failed for "${sampName}" at some point !!" && exit 1
    ################

    echo 'MergeVCF '${sampName}' files seems to be completed successfully' $(date)

    rm ${WD}/genotyping/${sampName}/${sampName}*_gt.log

  #### --------- Filter VCF --------- ####

    for q in ${quality[@]}
    do

      echo 'Running filters on '${sampName}' files with Q'${q} $(date)

      java -Xmx6g -jar ${NGSEP} FilterVCF -q ${q} ${WD}/genotyping/${sampName}/${sampName}.vcf 1> ${WD}/genotyping/${sampName}/${sampName}_q${q}.vcf

      ################
      test ! -s ${WD}/genotyping/${sampName}/${sampName}_q${q}.vcf && echo "Error: Filter failed for "${sampName}" at some point !!" && exit 1
      ################

      echo 'Filters on '${sampName}' VCF seems to be completed successfully' $(date)

  #### -------- Compare VCF --------- ####

      echo 'Comparing VCF for '${sampName} $(date)
  
      java -Xmx6g -jar ${NGSEP} CompareVCF -g 0 -d 100 ${referenceIDs} ${WD}/genotyping/${sampName}/${sampName}_q${q}.vcf ${WD}/genotyping/${sampName}/${sampName}_q${q}.vcf > ${WD}/genotyping/${sampName}/${sampName}_CompareVCF_q${q}.txt

      ################
      test ! -s ${WD}/genotyping/${sampName}/${sampName}_CompareVCF_q${q}.txt && echo "Error: CompareVCF failed for "${sampName}" at some point !!" && exit 1
      ################

    done

    rm ${WD}/genotyping/${sampName}/${sampName}.vcf ${WD}/genotyping/${sampName}/${sampName}_q*.vcf ${WD}/genotyping/${sampName}/${sampName}*_gt.vcf

  #### -------- Merge BAMs --------- ####

    echo 'Starting Merge BAMs on '${sampName}' files' $(date)

    ${samtoolsExec}/samtools merge ${WD}/genotyping/${sampName}/${sampName}_bowtie2_sorted.bam	${WD}/genotyping/${sampName}/${sampName}*_bowtie2_sorted.bam

    ################
    test ! -s ${WD}/genotyping/${sampName}/${sampName}_bowtie2_sorted.bam && echo "Error: Merge BAMs failed for "${sampName}" at some point !!" && exit 1
    ################

    echo 'Merge BAMs on '${sampName}' files seems to be completed successfully' $(date)

  #### --- FindVariants for combined BAM ----- ####

    echo 'Starting FindVariants from the combined BAMs for '${sampName}' sample' $(date)

    java -Xmx3g -jar ${NGSEP} FindVariants -h 0.0001 -maxBaseQS 30 -minQuality 40 -noRep -noRD -noRP -maxAlnsPerStartPos 100 -ignore5 ${i5} -ignore3 ${i3} -sampleId ${sampName} -knownSTRs ${STRs} ${REF} ${WD}/genotyping/${sampName}/${sampName}_bowtie2_sorted.bam ${WD}/genotyping/${sampName}/${sampName}_bowtie2_NGSEP >& ${WD}/genotyping/${sampName}/${sampName}_bowtie2_NGSEP.log
    bgzip ${WD}/genotyping/${sampName}/${sampName}_bowtie2_NGSEP.vcf

    ################
    if [[ ! `tail -1 ${WD}/genotyping/${sampName}/${sampName}_bowtie2_NGSEP.log` == *Completed* ]]; then echo "Error: Find variants failed at some point for "${sampName}; exit 1; fi 
    ################

  done ) &

done
wait

rm ${WD}/genotyping/*tmp

echo ''; echo 'CompareRepeatedSamples for '${popName}' files seems to be completed' $(date)