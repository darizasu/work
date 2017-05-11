#!/bin/bash

  # Key variables to specify

# This is the working directory full path. It should contain two directories: 
# 'reads' and 'mapping'. 'reads' must have a subdirectory called 'lane', 
# which contains the raw reads.
WD=/bioinfo2/projects/bean_tmp/JuanDavid/bodo_tests/GBSRRR

# This is your population's name
popName=RRR_introgressions

# The number of subprocesses you want to run. It depends on the number of available cores.
numThreads=18

# Specify the task(s) you want to perform. Include only the initial capital letter in a single string.
# It can include 'M'erge variants, 'F'ind variants, 'V'CF Merge, 'A'nnotate, 'R'un filters.
# To run all these analyses: ./runPopulation.sh 'MFVAR'
TASKS=$1

# This file contains the sample ID to be removed from the VCF file.
# Take a look at NGSEP FilterVCF -saf -fs flags
# This file should be located at ${WD}/population, otherwise its path must be specified.
samples2remove=/bioinfo2/projects/bean_tmp/JuanDavid/bodo_tests/GBSRRR/population/samples2remove_RRR.txt

# Specify the quality thresholds for FilterVCF
quality=(20 40 60)

# The file specified in the 'samples2merge' parameter contains info about 
# the location of BAM and VCF files for every sample.
# It also contains some of the parameters to call variants using NGSEP.
# Its full path must be specified. This file may contain comment lines starting with '#'
# This 'samples2merge' file contains 4 tab-separated columns:
# /path/to/sampleID[tab]sample_name[tab]ignore5[tab]ignore3

  # 1) '/path/to/sample' gives the full location + sample-prefix to a sample's VCF and BAM files
  # In other words, you must be able to find these files in the specified directory:
  #   /path/to/sample_bowtie2_NGSEP.vcf.gz
  #   /path/to/sample_bowtie2_sorted.bam
  # These files are produced after running 'runPlate.sh' correctly.

  # 2) 'sample_name' is the name that every sample will take in the final VCF file.
  # Be aware to avoid repeated names in this second column.
  # In case of repeated samples, you must name every sample as (e.g.): 
  # sample_name-p01F12, the '-p' is MANDATORY after 'samplename' ('p' stands for 'plate'),
  # and p01F12 may indicate the plate number and well that contained that sample.

  # 3-4) ignore5 and ignore3 are to ignore this many base pairs from the
  # 5' and 3' ends in NGSEP FindVariants.
 
# See the example below:
# /bioinfo2/projects/GBSplates/01/mapping/ALB_213	ALB_213-p01H10	4	10
samples2merge=/bioinfo2/projects/bean_tmp/JuanDavid/bodo_tests/GBSRRR/Samples_Abril2017.txt

# In case you DO NOT want to run MergeVariants with the whole list of samples
# specified in 'samples2merge', please use the parameter 'myVariants' to specify
# an empty VCF file with its full path. This file is produced after 
# running NGSEP MergeVariants with other set of VCF files.
myVariants=/bionas0/bean/populations/RRR/genotyping/Variants_RRR_noGenotype.vcf

  # Path to Software used

NGSEP=/data/software/NGSEPcore_3.0.2.jar

  # Reference genome files

# The reference genome file, in fasta format and indexed using bowtie2-index
REF=/data/references/bean/v2.1/bowtie2/Pvulgaris_442_v2.0.fa

# The annotation file for the reference genome provided previously, in gff3 format.
REFGFF=/data/references/bean/v2.1/annotation/Pvulgaris_v2_genes.gff3

# A file with a list of repetitive regions in the reference genome.
# Check NGSEP-FilterVCF for more details.
REPS=/data/references/bean/v2.1/repeats/Pvulgaris_v2_repMasked.list

# A file with the list of sequence identifiers from the reference genome, one by line.
refIDs=/data/references/bean/v2.1/assembly/Pvulgaris_442_v2.0_seqnames.txt


#### ---------------------------------------------------------------------- ####
#### ----------- DO NOT MODIFY ANYTHING FROM THIS POINT FORWARD ----------- ####
#### ---------------------------------------------------------------------- ####


echo -e '\nThis run was executed by:  '$(whoami)'\n'

# Get the 1st column from 'samples2merge' and put every line 
# as an element of a bash array called 'locations'. Avoid '#' lines
locations=(`grep -v '^#' ${samples2merge} | cut -f 1 | tr '\n' ' '`)

# Get the 2nd column from 'samples2merge' file and put every line 
# as an element of a bash array called 'list'. Avoid '#' lines
list=(`grep -v '^#' ${samples2merge} | cut -f 2 | tr '\n' ' '`)

echo -e 'This run contains the following samples:\n'${list[@]}'\n'

# The whole list of samples is divided in 'nThreads' sublists.
# Every sublist is run in the background. It is an individual file 
# called tmpList_XXX.tmp that contains chunks of the original list of samples.
# It won't continue until all the sublists have finished.
function assignThreads {
  for index in `seq ${#list[@]}`
    do samplesPerList=`expr ${index} % ${numThreads}`
    grep -v '#' ${samples2merge} | sed "${index}q;d" >> tmpList_${samplesPerList}.tmp
  done
}

# Create a symlink for specified files named as every element of ${list}
function createSymlynks {
  for index in ${!list[@]} 
  do ln -s ${locations[${index}]}_$1 ./${list[${index}]}_$1
  done
}

# Check the genotyping  and population paths exist
if [ ! -d ${WD}/genotyping ]
then mkdir ${WD}/genotyping
  echo 'The genotyping path '${WD}'/genotyping was created'
else if [ ! -d ${WD}/population ]
then mkdir ${WD}/population
  echo 'The population path '${WD}'/population was created'
fi

#### ------------------------------ ####
#### ------- Merge Variants ------- ####


if [[ ${TASKS} == *M* ]]
then

  cd ${WD}/genotyping

  echo -e '\nStarting merge variants on '${popName}' files\n' $(date)'\n'

  createSymlynks bowtie2_NGSEP.vcf.gz

  # VCF with the list of variants reported by the input files 
  # It doesn't contain any genotype information (Empty VCF)

  echo $(date) 'Merging variants from population'

  java -Xmx6g -jar ${NGSEP} MergeVariants ${refIDs} Variants_${popName}_noGenotype.vcf \
  ${WD}/genotyping/*.vcf.gz 2> ${popName}_mergevariants.log

  # Replace the content of the variable 'myVariants'
  myVariants=${WD}/genotyping/Variants_${popName}_noGenotype.vcf

  ################
  if [[ ! `tail -1 ${popName}_mergevariants.log` == *last* ]]
  then echo "Error: Merge variants failed at some point !!"; exit 1
  fi 
  ################

  rm ${WD}/genotyping/*.vcf.gz

  echo -e '\nMerge variants on '${popName}' files seems to be completed\n'$(date)'\n'

else

  echo -e '\nMerge Variants will not be executed this time\n'$(date)'\n'

fi


#### ------------------------------ ####
#### -------- Find Variants ------- ####


if [[ ${TASKS} == *F* ]]
then

  cd ${WD}/genotyping

  echo -e '\nStarting find variants on '${popName}' files\n'$(date)'\n'

  echo 'Total number of samples: '${#list[@]}

  assignThreads

  createSymlynks bowtie2_sorted.bam

  myNum=0

  for tmpFile in tmpList*tmp
  do

    # From the tmpList_XXX.tmp file (which is a chunk of samples2merge)
    # get the 2nd , 3rd and 4th columns and put them in a bash array
    myList=(`cut -f 2 ${tmpFile} | tr '\n' ' '`)
    i5=(`cut -f 3 ${tmpFile} | tr '\n' ' '`)
    i3=(`cut -f 4 ${tmpFile} | tr '\n' ' '`)

    echo -e 'File '${tmpFile}' contains the following samples:\n'${myList[@]}
    myNum=`expr ${myNum} + ${#myList[@]}`; echo -e 'No. of samples assigned: '${myNum}'\n'

  ( for p in ${myList[@]}
    do sleep 5

      echo $(date) 'Genotyping population variants in '${p}

      java -Xmx3g -jar ${NGSEP} FindVariants -h 0.0001 -maxBaseQS 30 -minQuality 0 \
      -noRep -noRD -noRP -maxAlnsPerStartPos 100 -ignore5 ${i5} -ignore3 ${i3} \
      -sampleId ${p} -knownVariants ${myVariants} ${REF} ./${p}_bowtie2_sorted.bam \
      ${p}_bowtie2_NGSEP_gt >& ${p}_bowtie2_NGSEP_gt.log

    ################
    if [[ ! `tail -1 ${p}_bowtie2_NGSEP_gt.log` == *Completed* ]]
    then echo "Error: FindVariants failed at some point for "${p}; exit 1
    fi
    ################

    done ) &

  done
  wait

  rm *tmp  *_bowtie2_sorted.ba*

  echo -e '\nGenotyping population variants on '${popName}' files seems to be completed\n'$(date)'\n'

else

  echo -e '\nFind Variants will not be executed this time\n'$(date)'\n'

fi


#### ------------------------------ ####
#### --------- Merge VCF ---------- ####


if [[ ${TASKS} == *V* ]]
then

  cd ${WD}/genotyping

  echo -e '\nMerging VCF files in '${popName}'\n'$(date)'\n'

  java -Xmx10g -jar ${NGSEP} MergeVCF ${refIDs} ./*_gt.vcf > ${WD}/population/${popName}.vcf

  ################
  test ! -s ${WD}/population/${popName}.vcf \
  && echo "Error: Merge failed at some point !!" && exit 1
  ################

  bgzip ${WD}/population/${popName}.vcf
  tabix -p vcf ${WD}/population/${popName}.vcf.gz

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
  && bgzip ${WD}/population/${popName}.vcf \
  && tabix -p vcf ${WD}/population/${popName}.vcf.gz

  java -Xmx6g -jar ${NGSEP} Annotate ./${popName}.vcf.gz \
  ${REFGFF} ${REF} > ${popName}_annotated.vcf

  ################
  test ! -s ${popName}_annotated.vcf \
  && echo "Error: Annotate failed at some point !!" && exit 1
  ################

  bgzip ${popName}_annotated.vcf
  tabix -p vcf ${popName}_annotated.vcf.gz

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

  echo -e '\nAnnotate will not be executed this time\n'$(date)'\n'

fi


#### ------------------------------ ####
#### --------- FilterVCF ---------- ####


if [[ ${TASKS} == *R* ]]
then

  cd ${WD}/population

  echo -e '\nRunning filters on  '${popName}'  files\n'$(date)'\n'

  for q in ${quality[@]}
  do

  ( java -Xmx6g -jar ${NGSEP} FilterVCF -saf ${samples2remove} \
    -fs -q ${q} ${popName}_annotated.vcf.gz 1> ${popName}_annotated_q${q}.vcf

    java -Xmx3g -jar ${NGSEP} SummaryStats -m 0 \
    ${popName}_annotated_q${q}.vcf > ${popName}_annotated_q${q}_summary.stats

    ################
    test ! -s ${popName}_annotated_q${q}_summary.stats \
    && echo 'Error: FilterVCF for '${popName}'_annotated_q'${q}'.vcf failed at some point !!' \
    && exit 1
    ################

    bgzip ${popName}_annotated_q${q}.vcf
    tabix -p vcf ${popName}_annotated_q${q}.vcf.gz

    java -Xmx6g -jar ${NGSEP} FilterVCF -s -frs ${REPS} -fi -minMAF 0.05 -maxOH 0.06 \
    ${popName}_annotated_q${q}.vcf.gz 1> ${popName}_annotated_repMasked_q${q}_s_fi_maf05_oh06.vcf

    java -Xmx3g -jar ${NGSEP} SummaryStats -m 1 ${popName}_annotated_repMasked_q${q}_s_fi_maf05_oh06.vcf \
    > ${popName}_annotated_repMasked_q${q}_s_fi_maf05_oh06_summary.stats

    ################
    test ! -s ${popName}_annotated_repMasked_q${q}_s_fi_maf05_oh06_summary.stats \
    && echo 'Error: FilterVCF for '${popName}'_annotated_q'${q}'_s_fi_maf05_oh06.vcf failed at some point !!' \
    && exit 1
    ################

    bgzip ${popName}_annotated_repMasked_q${q}_s_fi_maf05_oh06.vcf
    tabix -p vcf ${popName}_annotated_repMasked_q${q}_s_fi_maf05_oh06.vcf.gz

  ) &

  done
  wait

  echo -e '\nFilters on '${popName}' files seems to be completed\n'$(date)'\n'

else

echo -e '\nFilterVCF will not be executed this time\n'$(date)'\n'

fi
