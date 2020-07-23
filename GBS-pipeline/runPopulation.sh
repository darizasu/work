#!/bin/bash

  # Key variables to specify

# This is the working directory full path. It creates the directories 'genotyping' 
# and 'population' in case they don't exist.
WD=/bioinfo1/projects/bean/joint_pops/ADP_VEC_VEF/v2.1

# This is your population's name
popName=ADP_MIP_VEC_VEF

# The number of subprocesses you want to run. It depends on the number of available cores.
numThreads=20

# Specify the task(s) you want to perform. Include only the initial capital letter in a single string.
# It can include 'V'ariantsMerge, 'G'enotyping, 'M'ergeVCF, 'A'nnotate, 'F'ilterVCF.
# To run all these analyses: ./runPopulation.sh 'VGMAF'
TASKS=$1

# This file contains the sample ID to be removed from the VCF file.
# Take a look at NGSEP FilterVCF -saf -fs flags
# The full path to this file must be specified.
samples2remove=

# Specify the quality thresholds for FilterVCF
quality=(40)

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
samples2population=/bioinfo1/projects/bean/joint_pops/ADP_VEC_VEF/v2.1/samples_in_ADP_MIP_VEC_VEF.txt

# In case you DO NOT want to run MergeVariants with the whole list of samples
# specified in 'samples2population', please use the parameter 'myVariants' to specify
# a list of variants in VCF format with its full path. This file is produced after 
# running NGSEP MergeVariants with other set of VCF files.
myVariants=

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

NGSEP=/home/dariza/software/NGSEP/NGSEPcore_4.0.0.jar
bgzip=/usr/bin/bgzip
tabix=/usr/bin/tabix

  # Reference genome files

# The reference genome file, in fasta format indexed using bowtie2-index
REF=/bioinfo1/references/bean/G19833/v2.1/bowtie2/Pvulgaris_442_v2.0.fa

# The annotation file for the reference genome provided previously, in gff3 format.
REFGFF=/bioinfo1/references/bean/G19833/v2.1/annotation/Pvulgaris_v2_genes.gff3

# A file with a list of repetitive regions in the reference genome.
# Check NGSEP-FilterVCF for more details.
REPS=/bioinfo1/references/bean/G19833/v2.1/repeats/Pvulgaris442_repmasked.txt


#### ---------------------------------------------------------------------- ####
#### ----------- DO NOT MODIFY ANYTHING FROM THIS POINT FORWARD ----------- ####
#### ---------------------------------------------------------------------- ####


echo -e '\nThis run was executed by:  '$(whoami)'\n'

# Make sure you dont have ASCII text with CRLF line terminators
dos2unix ${samples2remove} ${samples2population}

# Create the refIDs file from the reference genome
grep '^>' ${REF} | sed "s/>//g" > ${WD}/Pvulgaris_seqnames.txt
refIDs=${WD}/Pvulgaris_seqnames.txt

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

  java -Xmx6g -jar ${NGSEP} MergeVariants -s ${refIDs} -o Variants_${popName}_noGenotype.vcf \
  ${WD}/genotyping/*${VCFext} 2> ${popName}_mergevariants.log

  # Replace the content of the variable 'myVariants'
  myVariants=${WD}/genotyping/Variants_${popName}_noGenotype.vcf

  ################
  if [[ ! `tail -1 ${popName}_mergevariants.log` == *Merged* ]]
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
    do sleep 2

      echo $(date) 'Genotyping population variants in '${myList[${index}]}

      java -Xmx20g -jar ${NGSEP} SingleSampleVariantsDetector -h 0.0001 -maxBaseQS 30 -minQuality 0 \
      -maxAlnsPerStartPos 100 -ignore5 ${Ifive[${index}]} -ignore3 ${Ithree[${index}]} \
      -sampleId ${myList[${index}]} -knownVariants ${myVariants} -r ${REF} -i ./${myList[${index}]}_${BAMext} \
      -o ${myList[${index}]}_${gtVCFext} >& ${myList[${index}]}_${gtVCFext}.log

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

  java -Xmx30g -jar ${NGSEP} MergeVCF -o /dev/stdout -s ${refIDs} ./*_${gtVCFext}.vcf.gz | \
  ${bgzip} > ${WD}/population/${popName}.vcf.gz

  ################
  test ! -s ${WD}/population/${popName}.vcf.gz \
  && echo "Error: Merge failed at some point !!" && exit 1
  ################

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

  java -Xmx6g -jar ${NGSEP} Annotate -i ./${popName}.vcf.gz \
  -t ${REFGFF} -r ${REF} -o /dev/stdout | ${bgzip} > ${popName}_annotated.vcf.gz

  ################
  test ! -s ${popName}_annotated.vcf.gz \
  && echo "Error: Annotate failed at some point !!" && exit 1
  ################

  ${tabix} -p vcf ${popName}_annotated.vcf.gz

  echo -e '\nVCF file in '${popName}' was annotated\n'$(date)'\n'

  # Summary Statistics
 
  echo $(date) 'Summary statistics for '${popName}'_annotated'

  java -Xmx6g -jar ${NGSEP} SummaryStats \
  -i ${popName}_annotated.vcf.gz -o ${popName}_annotated_summary.stats

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

  test -e ${WD}/population/${popName}_annotated.vcf.gz && ext='_annotated' || ext=''

  echo -e '\nRunning filters on '${popName}' files\n'$(date)'\n'

  for q in ${quality[@]}
  do

  ( java -Xmx6g -jar ${NGSEP} VCFFilter -saf ${samples2remove} \
    -fs -q ${q} -i ${popName}${ext}.vcf.gz -o /dev/stdout | ${bgzip} 1> ${popName}${ext}_q${q}.vcf.gz

    java -Xmx3g -jar ${NGSEP} SummaryStats \
    -i ${popName}${ext}_q${q}.vcf.gz -o ${popName}${ext}_q${q}_summary.stats

    ################
    test ! -s ${popName}${ext}_q${q}_summary.stats \
    && echo 'Error: FilterVCF for '${popName}${ext}'_q'${q}'.vcf.gz failed at some point !!' \
    && exit 1
    ################

    ${tabix} -p vcf ${popName}${ext}_q${q}.vcf.gz

    java -Xmx6g -jar ${NGSEP} VCFFilter -s -frs ${REPS} -fi -minMAF 0.05 -maxOH 0.05 \
    -i ${popName}${ext}_q${q}.vcf.gz -o /dev/stdout | ${bgzip} 1> ${popName}${ext}_repMasked_q${q}_s_fi_maf05_oh05.vcf.gz

    java -Xmx3g -jar ${NGSEP} SummaryStats -m 1 -i ${popName}${ext}_repMasked_q${q}_s_fi_maf05_oh05.vcf.gz \
    -o ${popName}${ext}_repMasked_q${q}_s_fi_maf05_oh05_summary.stats

    ################
    test ! -s ${popName}${ext}_repMasked_q${q}_s_fi_maf05_oh05_summary.stats \
    && echo 'Error: FilterVCF for '${popName}${ext}'_q'${q}'_s_fi_maf05_oh05.vcf.gz failed at some point !!' \
    && exit 1
    ################

    ${tabix} -p vcf ${popName}${ext}_repMasked_q${q}_s_fi_maf05_oh05.vcf.gz

  ) &

  done
  wait

  echo -e '\nFilters on '${popName}' files seems to be completed\n'$(date)'\n'

else

  echo -e '\nFilterVCF will not be executed this time\n'$(date)'\n'

fi

rm ${refIDs}
