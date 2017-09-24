#!/bin/bash

  # Key variables to specify

# This is the working directory full path. It creates the directories 'genotyping' 
# and 'population' in case they don't exist.
WD=/bioinfo1/projects/bean/all_genotypes

# This is your population's name
popName=all_bean_summary

# The number of subprocesses you want to run. It depends on the number of available cores.
numThreads=6

# Specify the task(s) you want to perform. Include only the initial capital letter in a single string.
# It can include 'V'ariantsMerge, 'G'enotyping, 'M'ergeVCF, 'A'nnotate, 'F'ilterVCF.
# To run all these analyses: ./runPopulation.sh 'VGMAF'
TASKS=$1

# This file contains the sample ID to be removed from the VCF file.
# Take a look at NGSEP FilterVCF -saf -fs flags
# The full path to this file must be specified.
samples2remove=/bioinfo1/projects/bean/all_genotypes/samples2remove_all_bean_summary.txt

# Specify the quality thresholds for FilterVCF
quality=(20 40 60)

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
  # Be aware to avoid repeated names in this second column.
  # In case of repeated samples, you must name every sample as (e.g.): 
  # sample_name-p01F12, the '-p' is mandatory (MANDATORY) after 'samplename' ('p' stands for 
  # 'plate'), and p01F12 may indicate the plate number and well that contained that sample.

  # 3-4) ignore5 and ignore3 are to ignore this many base pairs from the
  # 5' and 3' ends in NGSEP FindVariants.
 
# See the example below:
# /bioinfo2/projects/GBSplates/01/mapping/ALB_213	ALB_213-p01H10	4	10
samples2population=/bioinfo1/projects/bean/all_genotypes/samples_in_all_bean_summary.txt

# In case you DO NOT want to run MergeVariants with the whole list of samples
# specified in 'samples2population', please use the parameter 'myVariants' to specify
# a list of variants in VCF format with its full path. This file is produced after 
# running NGSEP MergeVariants with other set of VCF files.
myVariants=/bioinfo1/projects/bean/WGS_tpg/Variants_WGS_tpg_NO_phred33.vcf.gz

# The following variables define the file extension that the BAM and VCF files have in 
# the location specified in the '/path/to/sample' of 'samples2population' (1st column). 
# For example, if BAMext = bowtie2_sorted.bam and VCFext = bowtie2_NGSEP.vcf, then the 
# '/path/to/sample' of 'samples2population' should take you to 
# '/path/to/sample_bowtie2_sorted.bam' and '/path/to/sample_bowtie2_NGSEP.vcf'
BAMext=bowtie2_sorted.bam
VCFext=bowtie2_NGSEP.vcf.gz

# The following variable define the extension that the VCF file produced after 
# Genotyping will have. In other words, if gtVCFext = bowtie2_NGSEP_gt then
# you will find all genotyping files as ${WD}/genotyping/sample_bowtie2_NGSEP.vcf
gtVCFext=bowtie2_NGSEP_gt

  # Path to Software used

NGSEP=/data/software/NGSEPcore_3.0.2.jar
bgzip=/usr/bin/bgzip
tabix=/usr/bin/tabix

  # Reference genome files

# The reference genome file, in fasta format indexed using bowtie2-index
REF=/data/references/bean/v2.1/bowtie2/Pvulgaris_442_v2.0.fa

# The annotation file for the reference genome provided previously, in gff3 format.
REFGFF=/data/references/bean/v2.1/annotation/Pvulgaris_v2_genes.gff3

# A file with a list of repetitive regions in the reference genome.
# Check NGSEP-FilterVCF for more details.
REPS=/bionas1/bean/datasetPapers/201707_WGS_64vars/reference/Pvulgaris442_repmasked.txt

# A file with the list of sequence identifiers from the reference genome, one by line.
refIDs=/data/references/bean/v2.1/assembly/Pvulgaris_442_v2.0_seqnames.txt


#### ---------------------------------------------------------------------- ####
#### ----------- DO NOT MODIFY ANYTHING FROM THIS POINT FORWARD ----------- ####
#### ---------------------------------------------------------------------- ####


echo -e '\nThis run was executed by:  '$(whoami)'\n'

# Get the 1st column from 'samples2population' and put every line 
# as an element of a bash array called 'locations'. Avoid '#' lines
locations=(`grep -v '^#' ${samples2population} | cut -f 1 | tr '\n' ' '`)

# Get the 2nd column from 'samples2population' file and put every line 
# as an element of a bash array called 'list'. Avoid '#' lines
list=(`grep -v '^#' ${samples2population} | cut -f 2 | tr '\n' ' '`)

echo -e 'This run contains '${#list[@]}' samples:\n'${list[@]}'\n'

# The whole list of samples is divided in 'nThreads' sublists.
# Every sublist is run in the background. It is an individual file 
# called tmpList_XXX.tmp that contains chunks of the original list of samples.
# It won't continue until all the sublists have finished.
function assignThreads {
  for index in `seq ${#list[@]}`
    do samplesPerList=`expr ${index} % ${numThreads}`
    grep -v '#' ${samples2population} | sed "${index}q;d" >> tmpList_${samplesPerList}.tmp
  done
}

# Create a symlink for specified files named as every element of ${list}
function createSymlinks {
  for index in ${!list[@]} 
  do ln -s ${locations[${index}]}_$1 ./${list[${index}]}_$1
  done
}

# Check the genotyping  and population paths exist
if [ ! -d ${WD}/genotyping ]
then mkdir ${WD}/genotyping
  echo 'The genotyping path '${WD}'/genotyping was created'
fi
if [ ! -d ${WD}/population ]
then mkdir ${WD}/population
  echo 'The population path '${WD}'/population was created'
fi

#### ------------------------------ ####
#### ------- Merge Variants ------- ####


if [[ ${TASKS} == *V* ]]
then

  cd ${WD}/genotyping

  echo -e '\nStarting merge variants on '${popName}' files\n' $(date)'\n'

  createSymlinks ${VCFext}

  # VCF with the list of variants reported by the input files 
  # It doesn't contain any genotype information (Empty VCF)

  echo $(date) 'Merging variants from population'

  java -Xmx6g -jar ${NGSEP} MergeVariants ${refIDs} Variants_${popName}_noGenotype.vcf \
  ${WD}/genotyping/*${VCFext} 2> ${popName}_mergevariants.log

  # Replace the content of the variable 'myVariants'
  myVariants=${WD}/genotyping/Variants_${popName}_noGenotype.vcf

  ################
  if [[ ! `tail -1 ${popName}_mergevariants.log` == *last* ]]
  then echo "Error: Merge variants failed at some point !!"; exit 1
  fi 
  ################

  rm ${WD}/genotyping/*${VCFext} ${popName}_mergevariants.log

  echo -e '\nMerge variants on '${popName}' files seems to be completed\n'$(date)'\n'

else

  echo -e '\nMerge Variants will not be executed this time\n'$(date)'\n'

fi


#### ------------------------------ ####
#### -------- Find Variants ------- ####


if [[ ${TASKS} == *G* ]]
then

  cd ${WD}/genotyping

  echo -e '\nStarting find variants on '${popName}' files\n'$(date)'\n'

  echo 'Total number of samples: '${#list[@]}

  assignThreads

  createSymlinks ${BAMext}

  myNum=0

  for tmpFile in tmpList*tmp
  do

    # From the tmpList_XXX.tmp file (which is a chunk of samples2population)
    # get the 2nd , 3rd and 4th columns and put them in a bash array
    myList=(`cut -f 2 ${tmpFile} | tr '\n' ' '`)
    Ifive=(`cut -f 3 ${tmpFile} | tr '\n' ' '`)
    Ithree=(`cut -f 4 ${tmpFile} | tr '\n' ' '`)

    echo -e 'File '${tmpFile}' contains '${#myList[@]}' samples:\n'${myList[@]}
    myNum=`expr ${myNum} + ${#myList[@]}`; echo -e 'No. of samples assigned: '${myNum}'\n'

  ( for index in ${!myList[@]}
    do sleep 5

      echo $(date) 'Genotyping population variants in '${myList[${index}]}

      java -Xmx3g -jar ${NGSEP} FindVariants -h 0.0001 -maxBaseQS 30 -minQuality 0 \
      -noRep -noRD -noRP -maxAlnsPerStartPos 100 -ignore5 ${Ifive[${index}]} -ignore3 ${Ithree[${index}]} \
      -sampleId ${myList[${index}]} -knownVariants ${myVariants} ${REF} ./${myList[${index}]}_${BAMext} \
      ${myList[${index}]}_${gtVCFext} >& ${myList[${index}]}_${gtVCFext}.log

      ${bgzip} ${myList[${index}]}_${gtVCFext}.vcf

    done ) &

  done
  wait

  rm *tmp  *_${BAMext}

  ################
  # Check genotyping failures
  numErrors=0
  for index in ${!list[@]}
  do if [[ ! `tail -1 ${list[${index}]}_${gtVCFext}.log` == *Completed* ]]
    then numErrors=`expr ${numErrors} + 1`
    echo 'Error: FindVariants failed at some point for '${list[${index}]}
  fi; done
  if [[ ${numErrors} > 0 ]]; then echo 'Error: FindVariants failed for '${numErrors}' samples'; exit 1; fi
  ################

  echo -e '\nGenotyping population variants on '${popName}' files seems to be completed\n'$(date)'\n'

else

  echo -e '\nFind Variants will not be executed this time\n'$(date)'\n'

fi


#### ------------------------------ ####
#### --------- Merge VCF ---------- ####


if [[ ${TASKS} == *M* ]]
then

  cd ${WD}/genotyping

  echo -e '\nMerging VCF files in '${popName}'\n'$(date)'\n'

  java -Xmx10g -jar ${NGSEP} MergeVCF ${refIDs} ./*_${gtVCFext}.vcf.gz | \
  ${bgzip} > ${WD}/population/${popName}.vcf.gz

  ################
  test ! -s ${WD}/population/${popName}.vcf.gz \
  && echo "Error: Merge failed at some point !!" && exit 1
  ################

  # ${bgzip} ${WD}/population/${popName}.vcf
  ${tabix} -p vcf ${WD}/population/${popName}.vcf.gz

  echo -e '\nVCF files in '${popName}' were merged\n'$(date)'\n'

else

  echo -e '\nMerge will not be executed this time\n'$(date)'\n'

fi


#### ------------------------------ ####
#### -------- Annotate VCF -------- ####


if [[ ${TASKS} == *A* ]]
then

  cd ${WD}/population

  echo -e '\nAnnotating VCF file '${popName}'\n'$(date)'\n'

  test ! -s ${WD}/population/${popName}.vcf.gz \
  && ${bgzip} ${WD}/population/${popName}.vcf \
  && ${tabix} -p vcf ${WD}/population/${popName}.vcf.gz

  java -Xmx6g -jar ${NGSEP} Annotate ./${popName}.vcf.gz \
  ${REFGFF} ${REF} | ${bgzip} > ${popName}_annotated.vcf.gz

  ################
  test ! -s ${popName}_annotated.vcf.gz \
  && echo "Error: Annotate failed at some point !!" && exit 1
  ################

  # ${bgzip} ${popName}_annotated.vcf
  ${tabix} -p vcf ${popName}_annotated.vcf.gz

  echo -e '\nVCF file in '${popName}' was annotated\n'$(date)'\n'

  # Summary Statistics
 
  echo $(date) 'Summary statistics for '${popName}'_annotated'

  java -Xmx6g -jar ${NGSEP} SummaryStats -m 1 \
  ${popName}_annotated.vcf.gz > ${popName}_annotated_summary.stats

  ################
  test ! -s ${popName}_annotated_summary.stats \
  && echo "Error: SummaryStats failed at some point !!" && exit 1
  ################

else

  echo -e '\nAnnotation will not be executed this time\n'$(date)'\n'

fi


#### ------------------------------ ####
#### --------- FilterVCF ---------- ####


if [[ ${TASKS} == *F* ]]
then

  cd ${WD}/population

  echo -e '\nRunning filters on '${popName}' files\n'$(date)'\n'

  for q in ${quality[@]}
  do

  ( java -Xmx6g -jar ${NGSEP} FilterVCF -saf ${samples2remove} \
    -fs -q ${q} ${popName}_annotated.vcf.gz | ${bgzip} 1> ${popName}_annotated_q${q}.vcf.gz

    java -Xmx3g -jar ${NGSEP} SummaryStats -m 0 \
    ${popName}_annotated_q${q}.vcf.gz > ${popName}_annotated_q${q}_summary.stats

    ################
    test ! -s ${popName}_annotated_q${q}_summary.stats \
    && echo 'Error: FilterVCF for '${popName}'_annotated_q'${q}'.vcf.gz failed at some point !!' \
    && exit 1
    ################

    # ${bgzip} ${popName}_annotated_q${q}.vcf
    ${tabix} -p vcf ${popName}_annotated_q${q}.vcf.gz

    java -Xmx6g -jar ${NGSEP} FilterVCF -s -frs ${REPS} -fi -minMAF 0.05 -maxOH 0.06 \
    ${popName}_annotated_q${q}.vcf.gz | ${bgzip} 1> ${popName}_annotated_repMasked_q${q}_s_fi_maf05_oh06.vcf.gz

    java -Xmx3g -jar ${NGSEP} SummaryStats -m 1 ${popName}_annotated_repMasked_q${q}_s_fi_maf05_oh06.vcf.gz \
    > ${popName}_annotated_repMasked_q${q}_s_fi_maf05_oh06_summary.stats

    ################
    test ! -s ${popName}_annotated_repMasked_q${q}_s_fi_maf05_oh06_summary.stats \
    && echo 'Error: FilterVCF for '${popName}'_annotated_q'${q}'_s_fi_maf05_oh06.vcf.gz failed at some point !!' \
    && exit 1
    ################

    # ${bgzip} ${popName}_annotated_repMasked_q${q}_s_fi_maf05_oh06.vcf
    ${tabix} -p vcf ${popName}_annotated_repMasked_q${q}_s_fi_maf05_oh06.vcf.gz

  ) &

  done
  wait

  echo -e '\nFilters on '${popName}' files seems to be completed\n'$(date)'\n'

else

  echo -e '\nFilterVCF will not be executed this time\n'$(date)'\n'

fi
