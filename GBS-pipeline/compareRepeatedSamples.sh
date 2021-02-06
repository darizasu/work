#!/bin/bash

  # Key variables to specify

# This is the working directory full path. It creates the directories 'genotyping' 
# and 'population' in case they don't exist.
WD=/bioinfo1/projects/bean/duplicates

# This is your population's name
popName=duplicates

# The number of subprocesses you want to run. It depends on the number of available cores.
numThreads=16

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
samples2population=/bioinfo1/projects/bean/duplicates/samples_duplicates.txt


# The following variables define the file extension that the BAM and VCF files have in 
# the location specified in the '/path/to/sample' of 'samples2population' (1st column). 
# For example, if BAMext = bowtie2_sorted.bam and VCFext = bowtie2_NGSEP.vcf, then the 
# '/path/to/sample' of 'samples2population' should take you to 
# '/path/to/sample_bowtie2_sorted.bam' and '/path/to/sample_bowtie2_NGSEP.vcf'
BAMext=bowtie2_sorted.bam
VCFext=bowtie2_NGSEP
QCext=bowtie2_NGSEP_QC.vcf.gz

  # Path to Software used

NGSEP=/home/dariza/software/NGSEP/NGSEPcore_4.0.0.jar
samtools=/home/dariza/bin/samtools
bgzip=/usr/bin/bgzip
bcftools=/home/dariza/bin/bcftools

  # Reference genome files

REF=/bioinfo1/references/bean/Pvulgaris/G19833/v2.1/bowtie2/Pvulgaris_442_v2.0.fa
STRs=/bioinfo1/references/bean/Pvulgaris/G19833/v2.1/STRs/Pvulgaris_v2_strs.list


#### ---------------------------------------------------------------------- ####
#### ----------- DO NOT MODIFY ANYTHING FROM THIS POINT FORWARD ----------- ####
#### ---------------------------------------------------------------------- ####

echo -e '\nThis run was executed by:  '$(whoami)'\n'

# Make sure you dont have ASCII text with CRLF line terminators
dos2unix ${samples2population}

# Create the refIDs file from the reference genome
grep '^>' ${REF} | sed "s/>//g" > ${WD}/Pvulgaris_seqnames.txt
refIDs=${WD}/Pvulgaris_seqnames.txt

# Check the genotyping path exist

if [ ! -d ${WD}/genotyping ]
then mkdir ${WD}/genotyping
  echo -e 'The genotyping path '${WD}'/genotyping was created\n'
fi

cd ${WD}/genotyping

# Create a single file containing the samples to be included in every VCF.

# Take the whole list of samples
samples=(`grep -v '^#' ${samples2population} | cut -f 2 | tr '\n' ' '`)
# Get the list of duplicated samples from the whole list
dupSamples=(`grep -v '^#' ${samples2population} | cut -f 2 | sort -d | uniq -d | tr '\n' ' '`)
# Get the location for each sample's BAM and VCF files
locations=(`grep -v '^#' ${samples2population} | cut -f 1 | tr '\n' ' '`)

echo -e '\nThis run contains '${#dupSamples[@]}' samples sequenced more than once:\n' ${dupSamples[@]} '\n'

for index in ${!samples[@]}
do

# For each sample containing duplicates, put all its rows in an individual file named as the sample itself
  if [[ " ${dupSamples[@]} " =~ " ${samples[${index}]} " ]]
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


#### ------ Bring QC Variants ----- ####

cd ${WD}/genotyping

for file in ${WD}/genotyping/TMPsList*.tmp
do

  TMPs=(`cat ${file} | tr '\n' ' '`)
  echo -e 'File '${file}' contains '${#TMPs[@]}' samples:\n'${TMPs[@]}
  myNum=`expr ${myNum} + ${#TMPs[@]}`; echo -e 'No. of samples assigned: '${myNum}'\n'

 ( for sample in ${TMPs[@]}
  do  sleep 2

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
    sampName=`grep -v '^#' ${sample} | cut -f 2 | uniq`

    mkdir ${sampName}
    echo -e '\nThe sample '${sampName}' was sequenced '${#list[@]}' times:\n'${locations[@]}'\n' \
    >> ${WD}/genotyping/${sampName}/${sampName}_compRepeats.log

    for index in ${!list[@]}
    do

      ln -s ${locations[${index}]}*-p*${QCext} ./${sampName}/
      ln -s ${locations[${index}]}_${BAMext} ./${sampName}/${index}_${BAMext}

    done


  #### ---- Merge and Filter VCF ---- ####

    echo $(date) 'Merging and filtering VCF file for '${sampName} >> ${WD}/genotyping/${sampName}/${sampName}_compRepeats.log

    for QC in ${WD}/genotyping/${sampName}/*${QCext} ; do tabix -p vcf ${QC}; done

    ${bcftools} merge -Ou ${WD}/genotyping/${sampName}/*${QCext} | \
      ${bcftools} filter -Oz -o ${WD}/genotyping/${sampName}/${sampName}_q60.vcf.gz -S . -i 'GQ[*]>=60'

    ################
    if LC_ALL=C gzip -l ${WD}/genotyping/${sampName}/${sampName}_q60.vcf.gz | awk 'NR==2 {exit($2==0)}'
    then
      echo "Error: Merge or filter failed for "${sampName}" and Q 60 at some point !!" \
      >> ${WD}/genotyping/${sampName}/${sampName}_compRepeats.log
      continue
    fi
    ################


#### -------- Compare VCF --------- ####

    echo $(date) 'Comparing VCF file for '${sampName} >> ${WD}/genotyping/${sampName}/${sampName}_compRepeats.log

    java -Xmx6g -jar ${NGSEP} VCFComparator -g 0 -d 100 -r ${REF} \
    -i ${WD}/genotyping/${sampName}/${sampName}_q60.vcf.gz \
    -o ${WD}/genotyping/${sampName}/${sampName}_CompareVCF_q60.txt \
    >& ${WD}/genotyping/${sampName}/${sampName}_compRepeats.log

    ################
    test ! -s ${WD}/genotyping/${sampName}/${sampName}_CompareVCF_q60.txt \
    && echo "Error: CompareVCF failed for "${sampName}" at some point !!" \
    && echo "Error during execution of CompareVCF" >> ${WD}/genotyping/${sampName}/${sampName}_compRepeats.log \
    && continue
    ################

     rm ${WD}/genotyping/${sampName}/${sampName}*${QCext}*

  #### -------- Merge BAMs --------- ####

    echo $(date) 'Merging BAMs of '${sampName}' files' >> ${WD}/genotyping/${sampName}/${sampName}_compRepeats.log

    ${samtools} merge ${WD}/genotyping/${sampName}/${sampName}_${BAMext} \
    ${WD}/genotyping/${sampName}/*_${BAMext}

    ################
    test ! -s ${WD}/genotyping/${sampName}/${sampName}_${BAMext} \
    && echo "Error: Merge BAMs failed for "${sampName}" at some point !!" \
    && echo "Error during execution of Merge BAMs" >> ${WD}/genotyping/${sampName}/${sampName}_compRepeats.log \
    && continue
    ################


  #### --- FindVariants for combined BAM ----- ####

    echo $(date) 'Genotyping population variants from the combined BAMs for '${sampName} \
    >> ${WD}/genotyping/${sampName}/${sampName}_compRepeats.log

    java -Xmx3g -jar ${NGSEP} SingleSampleVariantsDetector \
    -h 0.0001 -maxBaseQS 30 -minQuality 40 -maxAlnsPerStartPos 100 \
    -ignore5 ${i5} -ignore3 ${i3} -sampleId ${sampName} -knownSTRs ${STRs} \
    -r ${REF} -i ${WD}/genotyping/${sampName}/${sampName}_${BAMext} \
    -o ${WD}/genotyping/${sampName}/${sampName}_${VCFext} \
    >& ${WD}/genotyping/${sampName}/${sampName}_${VCFext}.log
    ${bgzip} ${WD}/genotyping/${sampName}/${sampName}_${VCFext}.vcf

    ################
    if [[ ! `tail -1 ${WD}/genotyping/${sampName}/${sampName}_bowtie2_NGSEP.log` == *Completed* ]]
        then echo "Error: Find variants in "${sampName}" failed at some point !!" 
        echo "Error during execution of Find variants" >> ${WD}/genotyping/${sampName}/${sampName}_compRepeats.log
        continue
    fi 
    ################

    # if ! grep 'Error' ${WD}/genotyping/${sampName}/${sampName}_compRepeats.log
    #   then rm ${WD}/genotyping/${sampName}/${sampName}_compRepeats.log
    # fi

    echo -e ${WD}'/genotyping/'${sampName}/${sampName}'\t'${sampName}'\t'${i5}'\t'${i3} >> ${samples2population}

  done ) &

done
wait

rm ${WD}/genotyping/*tmp ${refIDs} ${WD}/genotyping/*/?_${BAMext} ${WD}/genotyping/*/*_bowtie2_NGSEP.log

echo -e '\nCompareRepeatedSamples for '${popName}' files seems to be completed successfully\n' $(date)