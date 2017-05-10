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

# The following file contains info about the location of BAM and VCF files for every sample.
# It also contains some of the parameters to call variants using NGSEP.
# Its full path must be specified.
# This file contains four tab-separated columns:
# /path/to/sample[tab]sample_name[tab]ignore5[tab]ignore3
# '/path/to/sample' gives the full location + sample prefix to a sample's VCF and BAM files
must contain a VCF and a BAM file for each sample on the same specified directory
# In this sense, you may be able to find the files /path/to/sample_bowtie2_NGSEP.vcf.gz and /path/to/sample_bowtie2_sorted.bam
# 'sample_name' is the name that the previously specified sample will take in the final VCF file. Be aware to avoid repeated names in this second column.
# In case of repeated samples, you must name every sample as (e.g.): sample_name-p01F12, the '-p' is MANDATORY after 'samplename' ('p' stands for 'plate'), and p01F12 may indicate the plate number and well that originated that sample.
# ignore5 is to ignore this many base pairs from the 5' end in NGSEP FindVariants
# ignore3 is to ignore this many base pairs from the 3' end in NGSEP FindVariants.
# See the example below:
# /bioinfo2/projects/GBSplates/01/mapping/ALB_213	ALB_213-p01H10	4	10
# This file may contain comment lines starting with '#'
samples2merge=/bioinfo2/projects/bean_tmp/JuanDavid/bodo_tests/GBSRRR/Samples_Abril2017.txt
myVariants=/bionas0/bean/populations/RRR/genotyping/Variants_RRR_noGenotype.vcf; # Specify an empty VCF file produced from NGSEP MergeVariants in case you do not want to run MergeVariants with the whole list of samples specified in 'samples2merge'. In other words, specify this option in case you don't run the 'F' option.

  # Path to Software used

NGSEP=/data/software/NGSEPcore_3.0.2.jar

  # Reference genome files

REF=/data/references/bean/v2.1/bowtie2/Pvulgaris_442_v2.0.fa
STRs=/data/references/bean/v2.1/strs/Pvulgaris_v2_strs.list
REFGFF=/data/references/bean/v2.1/annotation/Pvulgaris_v2_genes.gff3
REPS=/data/references/bean/v2.1/repeats/Pvulgaris_v2_repMasked.list
refIDs=/data/references/bean/v2.1/assembly/Pvulgaris_442_v2.0_seqnames.txt


#### ------------------------------------------------------------------------------------------------------------------------------- ####
#### ------------------------------------------------------------------------------------------------------------------------------- ####


locations=(`grep -v '^#' ${samples2merge} | cut -f 1 | tr '\n' ' '`); # Get the 1st column from samples2merge file and put every line as an element of a bash array called 'locations'. Avoid '#' lines
list=(`grep -v '^#' ${samples2merge} | cut -f 2 | tr '\n' ' '`); # Get the 2nd column from samples2merge file and put every line as an element of a bash array called 'locations'. Avoid '#' lines
echo 'This run was executed by:  '$(whoami)
echo 'This run contains the following samples:'; echo ${list[@]}


#### ------------------------------ ####
#### ------- Merge Variants ------- ####


if [[ ${TASKS} == *M* ]]
then

  if [ ! -d ${WD}/genotyping ]; then echo 'The genotyping path '${WD}'/mapping does not exist !'; mkdir ${WD}/genotyping; fi
  cd ${WD}/genotyping

  echo ''; echo 'Starting merge variants on '${popName}' files'; date; echo '';

 # Put all population VCFs in a single directory.

  for index in ${!list[@]}; do ln -s ${locations[${index}]}_bowtie2_NGSEP.vcf.gz ./${list[${index}]}_bowtie2_NGSEP.vcf.gz;  done; # Create a symlink for VCF files named as every element of ${locations}

 # VCF with the union of variants reported by the input files but without any genotype information (Empty VCF)

  echo $(date) 'Merging variants from population'

  java -Xmx6g -jar ${NGSEP} MergeVariants ${refIDs} Variants_${popName}_noGenotype.vcf ${WD}/genotyping/*.vcf.gz 2> ${popName}_mergevariants.log

  myVariants=${WD}/genotyping/Variants_${popName}_noGenotype.vcf

  ################
  if [[ ! `tail -1 ${popName}_mergevariants.log` == *last* ]]; then echo "Error: Merge variants failed at some point !!"; exit 1; fi 
  ################

  echo ''; echo 'Merge variants on '${popName}' files seems to be completed'; date; echo '';

else

  echo ''; echo 'Merge Variants will not be executed this time'; date; echo ''

fi


#### ------------------------------ ####
#### -------- Find Variants ------- ####


cd ${WD}/genotyping

if [[ ${TASKS} == *F* ]]
then

  echo ''; echo 'Starting find variants on  '${popName}'  files'; date; echo '';

  echo 'Total number of samples: '${#list[@]}

  # The whole list of samples is divided in 'nThreads' sublists,
  # and every sublist is run in the background. Every sublist 
  # is an individual file called tmpList_XXX.tmp
  # that contains chunks of the original 'samples2merge'.
  # It won't continue until all the sublists have finished.

  for index in `seq ${#list[@]}`; do samplesPerList=`expr ${index} % ${numThreads}`; grep -v '#' ${samples2merge} | sed "${index}q;d" >> tmpList_${samplesPerList}.tmp;  done

  for index in ${!list[@]}; do ln -s ${locations[${index}]}_bowtie2_sorted.bam ./${list[${index}]}_bowtie2_sorted.bam;  done; # Create a symlink for BAM files named as every element of 'locations'

  myNum=0

  for tmpFile in tmpList*tmp
  do

    myList=(`cut -f 2 ${tmpFile} | tr '\n' ' '`); # From the tmpList_XXX.tmp (which is a chunk of samples2merge) file get the 2nd column and put it in a bash array.
    i5=(`cut -f 3 ${tmpFile} | tr '\n' ' '`); # Get the 3rd column and put it in a bash array.
    i3=(`cut -f 4 ${tmpFile} | tr '\n' ' '`); # Get the 4th column and put it in a bash array.

    echo 'File '${tmpFile}' contains the following samples: '; echo ${myList[@]}

  ( for p in ${myList[@]}
    do

      sleep 5

      echo $(date) 'Genotyping population variants in '${p}

      java -Xmx3g -jar ${NGSEP} FindVariants -h 0.0001 -maxBaseQS 30 -minQuality 0 -noRep -noRD -noRP -maxAlnsPerStartPos 100 -ignore5 ${i5} -ignore3 ${i3} -sampleId ${p} -knownVariants ${myVariants} ${REF} ./${p}_bowtie2_sorted.bam ${p}_bowtie2_NGSEP_gt >& ${p}_bowtie2_NGSEP_gt.log

    ################
    if [[ ! `tail -1 ${p}_bowtie2_NGSEP_gt.log` == *Completed* ]]; then echo "Error: Find variants failed at some point for "${p}; exit 1; fi 
    ################

    done ) &

    myNum=`expr ${myNum} + ${#myList[@]}`; echo 'No. of samples assigned: '${myNum}; echo ''

  done
  wait

  rm *tmp  *_bowtie2_sorted.ba*

  echo ''; echo 'Genotyping population variants on '${popName}' files seems to be completed'; date; echo '';

else

  echo ''; echo 'Find Variants will not be executed this time'; date; echo ''

fi


#### ------------------------------ ####
#### --------- Merge VCF ---------- ####


cd ${WD}/genotyping

if [[ ${TASKS} == *V* ]]
then

  echo ''; echo 'Merging VCF files in '${popName}; date; echo '';

  java -Xmx10g -jar ${NGSEP} MergeVCF ${refIDs} ./*_gt.vcf > ${WD}/population/${popName}.vcf

  ################
  test ! -s ${WD}/population/${popName}.vcf && echo "Error: Merge failed at some point !!" && exit 1
  ################

  bgzip ${WD}/population/${popName}.vcf
  tabix -p vcf ${WD}/population/${popName}.vcf.gz

  echo ''; echo 'VCF files in '${popName}' were merged'; date; echo '';

else

  echo ''; echo 'Merge will not be executed this time'; date; echo ''

fi


#### ------------------------------ ####
#### -------- Annotate VCF -------- ####


cd ${WD}/population

if [[ ${TASKS} == *A* ]]
then

  echo ''; echo 'Annotating VCF file '${popName}; date; echo '';

  test ! -s ${WD}/population/${popName}.vcf.gz && bgzip ${WD}/population/${popName}.vcf && tabix -p vcf ${WD}/population/${popName}.vcf.gz

  java -Xmx6g -jar ${NGSEP} Annotate ./${popName}.vcf.gz ${REFGFF} ${REF} > ${popName}_annotated.vcf

  ################
  test ! -s ${popName}_annotated.vcf && echo "Error: Annotate failed at some point !!" && exit 1
  ################

  bgzip ${popName}_annotated.vcf
  tabix -p vcf ${popName}_annotated.vcf.gz

  echo ''; echo 'VCF file in '${popName}' was annotated'; date; echo '';

  # Summary Statistics
 
  echo $(date) 'Summary statistics for '${popName}'_annotated'

  java -Xmx6g -jar ${NGSEP} SummaryStats -m 1 ${popName}_annotated.vcf.gz > ${popName}_annotated_summary.stats

  ################
  test ! -s ${popName}_annotated_summary.stats && echo "Error: SummaryStats failed at some point !!" && exit 1
  ################

else

  echo ''; echo 'Annotate will not be executed this time'; date; echo ''

fi


#### ------------------------------ ####
#### --------- FilterVCF ---------- ####


echo ''; echo 'Running filters on  '${popName}'  files'; date; echo '';

cd ${WD}/population

if [[ ${TASKS} == *R* ]]
then

  for q in ${quality[@]}
  do

  ( sleep 5
    java -Xmx6g -jar ${NGSEP} FilterVCF -saf ${samples2remove} -fs -q ${q} ${popName}_annotated.vcf.gz 1> ${popName}_annotated_q${q}.vcf
    java -Xmx3g -jar ${NGSEP} SummaryStats -m 0 ${popName}_annotated_q${q}.vcf > ${popName}_annotated_q${q}_summary.stats

    ################
    test ! -s ${popName}_annotated_q${q}_summary.stats && echo 'Error: FilterVCF for '${popName}'_annotated_q'${q}'.vcf failed at some point !!' && exit 1
    ################

    bgzip ${popName}_annotated_q${q}.vcf
    tabix -p vcf ${popName}_annotated_q${q}.vcf.gz

    java -Xmx6g -jar ${NGSEP} FilterVCF -s -frs ${REPS} -minMAF 0.05 -maxOH 0.06 ${popName}_annotated_q${q}.vcf.gz 1> ${popName}_annotated_repMasked_q${q}_s_maf05_oh06.vcf
    java -Xmx3g -jar ${NGSEP} SummaryStats -m 1 ${popName}_annotated_repMasked_q${q}_s_maf05_oh06.vcf > ${popName}_annotated_repMasked_q${q}_s_maf05_oh06_summary.stats

    ################
    test ! -s ${popName}_annotated_repMasked_q${q}_s_maf05_oh06_summary.stats && echo 'Error: FilterVCF for '${popName}'_annotated_q'${q}'_s_maf05_oh06.vcf failed at some point !!' && exit 1
    ################

    bgzip ${popName}_annotated_repMasked_q${q}_s_maf05_oh06.vcf
    tabix -p vcf ${popName}_annotated_repMasked_q${q}_s_maf05_oh06.vcf.gz

  ) &

  done
  wait

  echo ''; echo 'Filters on '${popName}' files seems to be completed'; date; echo '';

else

  echo ''; echo 'FilterVCF will not be executed this time'; date; echo ''

fi
