#!/bin/bash

  # Key variables to specify

# This is the working directory full path. It creates the directories 'genotyping' 
# and 'population' in case they don't exist.
WD=/bioinfo1/projects/bean/VEF/v2.1

# This is your population's name
popName=VEF

# Specify the quality thresholds for FilterVCF
quality=(60)

# The number of subprocesses you want to run. It depends on the number of available cores.
numThreads=6

# The file specified in the 'samples2population' parameter contains info about 
# the location of BAM and VCF files for every sample.
# It also contains some of the parameters to call variants using NGSEP.
# Its full path must be specified. This file may contain comment lines starting with '#'
# This 'samples2population' file contains 4 tab-separated columns:
# /path/to/sampleID[tab]sample_name[tab]ignore5[tab]ignore3

  # 1) '/path/to/sample' gives the full location + sample-prefix to a sample's VCF and BAM files
  # In other words, you must be able to find these files in the specified directory 
  # (check parameters 'BAMext' and 'VCFext' below for more info):
  #   /path/to/sample_${VCFext}
  #   /path/to/sample_${BAMext}
  # These files are produced after running 'runPlate.sh' correctly.

  # 2) 'sample_name' is the name that every sample will take in the final VCF file.
  # WARNING ! Be aware to avoid repeated names in this second column.
  # In case of repeated samples, you must name every sample as (e.g.): 
  # sample_name-p01F12, the '-p' is mandatory (MANDATORY) after 'samplename' ('p' stands for 
  # 'plate'), and p01F12 may indicate the plate number and well that contained that sample.

  # 3-4) ignore5 and ignore3 are to ignore this many base pairs from the
  # 5' and 3' ends in NGSEP FindVariants.
 
# See the example below:
# /bioinfo2/projects/GBSplates/01/mapping/ALB_213	ALB_213-p01H10	4	10
samples2population=/bioinfo1/projects/bean/VEF/v2.1/samples_in_VEF.txt

# A file with the list of sequence identifiers from the reference genome, one by line.
referenceIDs=/data/references/bean/v2.1/assembly/Pvulgaris_442_v2.0_seqnames.txt

# Specify a list of predefined variants in a VCF file with its full path.
# This file is produced after running NGSEP MergeVariants with a 
# (preferably large) set of VCF files.
myVariants=/bioinfo1/projects/bean/WGS_tpg/Variants_WGS_tpg_NO_phred33.vcf

# The following variables define the file extension that the BAM and VCF files have in 
# the location specified in the '/path/to/sample' of 'samples2population' (1st column). 
# For example, if BAMext = bowtie2_sorted.bam and VCFext = bowtie2_NGSEP.vcf, then the 
# '/path/to/sample' of 'samples2population' should take you to 
# '/path/to/sample_bowtie2_sorted.bam' and '/path/to/sample_bowtie2_NGSEP.vcf'
BAMext=bowtie2_sorted.bam
VCFext=bowtie2_NGSEP

  # Path to Software used

NGSEP=/data/software/NGSEPcore_3.0.2.jar
samtoolsExec=/data/software/samtools/samtools-1.2
bgzip=/usr/bin/bgzip

  # Reference genome files

REF=/data/references/bean/v2.1/bowtie2/Pvulgaris_442_v2.0.fa
STRs=/data/references/bean/v2.1/strs/Pvulgaris_v2_strs.list


#### ---------------------------------------------------------------------- ####
#### ----------- DO NOT MODIFY ANYTHING FROM THIS POINT FORWARD ----------- ####
#### ---------------------------------------------------------------------- ####

echo -e '\nThis run was executed by:  '$(whoami)'\n'

# Check the genotyping path exist

if [ ! -d ${WD}/genotyping ]
then mkdir ${WD}/genotyping
  echo -e 'The genotyping path '${WD}'/genotyping was created\n'
fi

# Create a single file containing the samples to be included in every VCF.

cd ${WD}/genotyping

# Take the whole list of samples, and strip the -p*
samples=(`grep -v '^#' ${samples2population} | cut -f 2 | sed 's/-p.*//' | tr '\n' ' '`)
# Get the list of duplicated samples from the whole list
sampleIDs=(`grep -v '^#' ${samples2population} | cut -f 2 | sed 's/-p.*//' |  sort -d | uniq -d | tr '\n' ' '`)
# Get the location for each sample's BAM and VCF files
locations=(`grep -v '^#' ${samples2population} | cut -f 1 | tr '\n' ' '`)

echo -e '\nThis run contains '${#sampleIDs[@]}' samples sequenced more than once:\n' ${sampleIDs[@]} '\n'

for index in ${!samples[@]}
do

# For each sample containing duplicates, put all of its corresponding rows in an individual file named as the sample itself
  if [[ " ${sampleIDs[@]} " =~ " ${samples[${index}]} " ]]
  then
    index2=`expr ${index} + 1`
    grep -v '^#' ${samples2population} | sed -n -e ${index2}p >> ${samples[${index}]}.tmp
  else
    :
  fi

done

# Create an array of files for the tmp files created for each sample containing duplicates
myTMPs=(`ls *tmp`)

# Put every sample's tmp file in sublists according to the number of threads provided
for index in ${!myTMPs[@]}
do
  samplesPerList=`expr ${index} % ${numThreads}`
  echo ${myTMPs[${index}]} >> TMPsList_${samplesPerList}.tmp
done


#### -------- Find Variants ------- ####

cd ${WD}/genotyping

for file in ${WD}/genotyping/TMPsList*.tmp
do

  TMPs=(`cat ${file} | tr '\n' ' '`)
  echo -e 'File '${file}' contains '${#TMPs[@]}' samples:\n'${TMPs[@]}
  myNum=`expr ${myNum} + ${#TMPs[@]}`; echo -e 'No. of samples assigned: '${myNum}'\n'

 ( for sample in ${TMPs[@]}
  do  sleep 5

    echo -e $(date) 'Running all the analyses on '${sample}

    locations=(`grep -v '^#' ${sample} | cut -f 1 | tr '\n' ' '`)
    list=(`grep -v '^#' ${sample} | cut -f 2 | tr '\n' ' '`)
    # Put all ignore5 in a single array
    Ifive=(`grep -v '^#' ${sample} | cut -f 3 | tr '\n' ' '`)
    # Get the maximum ignore5 from the previous array
    i5=`echo ${Ifive[@]} | tr ' ' '\n' | sort -nr | head -1`
    # Put all ignore3 in a single array
    Ithree=(`grep -v '^#' ${sample} | cut -f 4 | tr '\n' ' '`)
    # Get the maximum ignore3 from the previous array
    i3=`echo ${Ithree[@]} | tr ' ' '\n' | sort -nr | head -1`
    sampName=`grep -v '^#' ${sample} | cut -f 2 | sed 's/-p.*//' | uniq`

    mkdir ${sampName}
    echo -e '\nThe sample '${sampName}' was sequenced '${#list[@]}' times:\n'${list[@]}'\n' \
    >> ${WD}/genotyping/${sampName}/${sampName}_compRepeats.log

    for index in ${!list[@]}
    do

      ln -s ${locations[${index}]}_${BAMext} ./${sampName}/${list[${index}]}_${BAMext}

      echo $(date) 'Genotyping population variants in '${list[${index}]} >> ${WD}/genotyping/${sampName}/${sampName}_compRepeats.log

      java -Xmx20g -jar ${NGSEP} FindVariants -h 0.0001 -maxBaseQS 30 \
      -minQuality 0 -noRep -noRD -noRP -maxAlnsPerStartPos 100 \
      -ignore5 ${Ifive[${index}]} -ignore3 ${Ithree[${index}]} \
      -sampleId ${list[${index}]} -knownVariants ${myVariants} ${REF} \
      ./${sampName}/${list[${index}]}_${BAMext} \
      ./${sampName}/${list[${index}]}_${VCFext}_gt \
      >& ./${sampName}/${list[${index}]}_${VCFext}_gt.log

      ################
      if [[ ! `tail -1 ./${sampName}/${list[${index}]}_${VCFext}_gt.log` == *Completed* ]]
        then echo "Error: Find variants in "${list[${index}]}" failed at some point !!" 
        echo "Error during execution of Find variants" >> ${WD}/genotyping/${sampName}/${sampName}_compRepeats.log
      fi 
      ################

    done


  #### --------- Merge VCF ---------- ####

    echo $(date) 'Merging VCF file for '${sampName} >> ${WD}/genotyping/${sampName}/${sampName}_compRepeats.log

    java -Xmx10g -jar ${NGSEP} MergeVCF ${referenceIDs} \
    ${WD}/genotyping/${sampName}/${sampName}*_gt.vcf \
    > ${WD}/genotyping/${sampName}/${sampName}.vcf

    ################
    test ! -s ${WD}/genotyping/${sampName}/${sampName}.vcf \
    && echo "Error: MergeVCF failed for "${sampName}" at some point !!" \
    && echo "Error during execution of MergeVCF" >> ${WD}/genotyping/${sampName}/${sampName}_compRepeats.log \
    && continue
    ################

    rm ${WD}/genotyping/${sampName}/${sampName}*_gt.log

  #### --------- Filter VCF --------- ####

    for q in ${quality[@]}
    do

      echo $(date) 'Running filters on '${sampName}' files with Q'${q} >> ${WD}/genotyping/${sampName}/${sampName}_compRepeats.log

      java -Xmx6g -jar ${NGSEP} FilterVCF -q ${q} \
      ${WD}/genotyping/${sampName}/${sampName}.vcf \
      1> ${WD}/genotyping/${sampName}/${sampName}_q${q}.vcf

      ################
      test ! -s ${WD}/genotyping/${sampName}/${sampName}_q${q}.vcf \
      && echo "Error: Filter failed for "${sampName}" at some point !!" \
      && echo "Error during execution of Filters with Q"${q} \
      >> ${WD}/genotyping/${sampName}/${sampName}_compRepeats.log \
      && continue
      ################


  #### -------- Compare VCF --------- ####

      echo $(date) 'Comparing VCF file for '${sampName} >> ${WD}/genotyping/${sampName}/${sampName}_compRepeats.log
  
      java -Xmx6g -jar ${NGSEP} CompareVCF -g 0 -d 100 ${referenceIDs} \
      ${WD}/genotyping/${sampName}/${sampName}_q${q}.vcf \
      ${WD}/genotyping/${sampName}/${sampName}_q${q}.vcf \
      > ${WD}/genotyping/${sampName}/${sampName}_CompareVCF_q${q}.txt

      ################
      test ! -s ${WD}/genotyping/${sampName}/${sampName}_CompareVCF_q${q}.txt \
      && echo "Error: CompareVCF failed for "${sampName}" at some point !!" \
      && echo "Error during execution of CompareVCF" >> ${WD}/genotyping/${sampName}/${sampName}_compRepeats.log \
      && continue
      ################

    done

    rm ${WD}/genotyping/${sampName}/${sampName}.vcf \
    ${WD}/genotyping/${sampName}/${sampName}_q*.vcf \
    ${WD}/genotyping/${sampName}/${sampName}*_gt.vcf

  #### -------- Merge BAMs --------- ####

    echo $(date) 'Merging BAMs of '${sampName}' files' >> ${WD}/genotyping/${sampName}/${sampName}_compRepeats.log

    ${samtoolsExec}/samtools merge ${WD}/genotyping/${sampName}/${sampName}_${BAMext} \
    ${WD}/genotyping/${sampName}/${sampName}*_${BAMext}

    ################
    test ! -s ${WD}/genotyping/${sampName}/${sampName}_${BAMext} \
    && echo "Error: Merge BAMs failed for "${sampName}" at some point !!" \
    && echo "Error during execution of Merge BAMs" >> ${WD}/genotyping/${sampName}/${sampName}_compRepeats.log \
    && continue
    ################


  #### --- FindVariants for combined BAM ----- ####

    echo $(date) 'Genotyping population variants from the combined BAMs for '${sampName} \
    >> ${WD}/genotyping/${sampName}/${sampName}_compRepeats.log

    java -Xmx3g -jar ${NGSEP} FindVariants -h 0.0001 -maxBaseQS 30 \
    -minQuality 40 -noRep -noRD -noRP -maxAlnsPerStartPos 100 \
    -ignore5 ${i5} -ignore3 ${i3} -sampleId ${sampName} -knownSTRs ${STRs} \
    ${REF} ${WD}/genotyping/${sampName}/${sampName}_${BAMext} \
    ${WD}/genotyping/${sampName}/${sampName}_${VCFext} \
    >& ${WD}/genotyping/${sampName}/${sampName}_${VCFext}.log
    ${bgzip} ${WD}/genotyping/${sampName}/${sampName}_${VCFext}.vcf

    ################
    if [[ ! `tail -1 ${WD}/genotyping/${sampName}/${sampName}_bowtie2_NGSEP.log` == *Completed* ]]
        then echo "Error: Find variants in "${sampName}" failed at some point !!" 
        echo "Error during execution of Find variants" >> ${WD}/genotyping/${sampName}/${sampName}_compRepeats.log
        continue
    fi 
    ################

    if ! grep 'Error' ${WD}/genotyping/${sampName}/${sampName}_compRepeats.log
      then rm ${WD}/genotyping/${sampName}/${sampName}_compRepeats.log
    fi

    echo -e ${WD}'/genotyping/'${sampName}/${sampName}'\t'${sampName}'\t'${i5}'\t'${i3} >> ${samples2population}

  done ) &

done
wait

rm ${WD}/genotyping/*tmp

echo -e '\nCompareRepeatedSamples for '${popName}' files seems to be completed successfully\n' $(date)
