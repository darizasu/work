#!/bin/bash

  # Key variables to specify

# This is the working directory full path. It should contain two directories: 
# 'reads' and 'mapping'. 'reads' must have a subdirectory called 'lane', 
# which contains the raw reads.
WD=/bioinfo2/projects/GBSplates/13

# This file should be located at ${WD}/reads/lane , otherwise its path must be specified. 
# Check NGSEP Deconvolute <INDEX_FILE> parameter for more info.
INDEXFILE=/bioinfo2/projects/GBSplates/13/reads/lane/barcodeMap_plate13.txt

# This file should be located at ${WD}/reads/lane , otherwise its path must be specified. 
# Check NGSEP Deconvolute -d flag for more info.
FILES2DECONV=/bioinfo2/projects/GBSplates/13/reads/lane/lanes_plate13.txt

# This is your plate's name
runName=plate_13

# The number of subprocesses you want to run. It depends on the number of available cores.
numThreads=20

# Specify the task(s) you want to perform. Include only the initial capital letter in a single string.
#It can include 'D'econvolution, 'T'rimming, 'M'apping, 'V'ariant-Discovery.
# To run all these analyses: ./runPlate.sh 'DTMV'
TASKS=$1

# Ignore this many base pairs from the 5' and 3' ends in NGSEP FindVariants.
# To get this numbers, run this script with the tasks 'DTM', then run the 'calculateReadPosGBS.sh' script 
# to get sequencing error bias for the entire plate and plot those results. Then decide the i5 
# and i3 parameters and specify them in the following lines. 
# Then run again this script with the task 'V'.
i5=7
i3=12

# This file must (MUST) be located at ${WD}/reads. This is a fasta file containing adapter 
# sequences to be removed from the deconvoluted reads. Check Trimmomatic manual for more info.
adapters=adapters_13-15.fa

  # Path to Software used

NGSEP=/data/software/NGSEPcore_3.0.2.jar
BOWTIE2=/data/software/bowtie2-2.3.0/bowtie2
PICARD=/data/software/picard-tools-1.140/picard.jar
Trimmomatic=/bioinfo1/software/Trimmomatic-0.36/trimmomatic-0.36.jar

  # Reference genome files

REF=/data/references/bean/v2.1/bowtie2/Pvulgaris_442_v2.0.fa
STRs=/data/references/bean/v2.1/strs/Pvulgaris_442_v2.0.fa.2.7.7.80.10.20.50.ngs


#### ---------------------------------------------------------------------- ####
#### ----------- DO NOT MODIFY ANYTHING FROM THIS POINT FORWARD ----------- ####
#### ---------------------------------------------------------------------- ####

echo -e '\nThis run was executed by:  '$(whoami)'\n'

if [ ! -d ${WD}/reads/lane ]; then echo 'The raw reads path '${WD}'/reads/lane does not exist !'; exit 1; fi
cd ${WD}/reads/lane

# Get sequencing run type, either single-end or paired-end from 
# the 'FILES2DECONV' file based on the number of columns it contains.
SEorPE=`head -1 ${FILES2DECONV} | tr '\t' '\n' | wc -l`

test ${SEorPE} == 4 && echo -e 'This run is working with Paired-end sequencing\n' || echo -e 'This run is working with Single-end sequencing\n'

# Get the list of samples from the 4th column in 'INDEXFILE' file, avoid first line. 
list=(`cut -f 4 ${INDEXFILE} | tail -n +2 | sort -u | tr '\n' ' '`)

echo -e 'This run contains the following samples:\n'${list[@]}'\n'

# The whole list of samples is divided in 'nThreads' sublists.
# Every sublist is run in the background. It is an individual file 
# called tmpList_XXX.tmp that contains chunks of the original list of samples.
# It won't continue until all the sublists have finished.
function assignThreads {
  for index in ${!list[@]}
  do  samplesPerList=`expr ${index} % ${numThreads}`
    echo ${list[${index}]} >> tmpList_${samplesPerList}.tmp
  done
}

# Check there are reads for every sample, exit otherwise.
function doWeHaveReads {
  numErrors=0
  for index in ${!list[@]}
  do if ! ls ${WD}/reads/${list[${index}]}*.fastq.gz > /dev/null 2>&1
    then numErrors=`expr ${numErrors} + 1`; echo 'Error: The sample '${list[${index}]}' does not have reads'
  fi; done
  if [[ ${numErrors} > 0 ]]; then echo 'Error: There are '${numErrors}' samples with no reads at '${WD}'/reads/'; exit 1; fi
}


#### ------------------------------ ####
#### ------- Deconvolution -------- ####

if [ ! -d ${WD}/reads ]; then echo 'The raw-reads path '${WD}'/reads/lane does not exist !'; exit 1; fi
cd ${WD}/reads/lane

if [[ ${TASKS} == *D* ]]
then

  echo -e '\nStarting deconvolution on '${runName}' files\n'$(date)'\n'

  java -jar ${NGSEP} Deconvolute -o ${WD}/reads -d ${FILES2DECONV} ${INDEXFILE}

fi

################
doWeHaveReads
################

echo -e '\nDeconvolution on '${runName}' files is done\n'$(date)'\n'


#### ------------------------------ ####
#### ----------- Trimming --------- ####


if [[ ${TASKS} == *T* ]]
then

  cd ${WD}/reads

  echo -e '\nStarting trimming on '${runName}' files'$(date)'\n'
  mkdir unTrimmed_reads
  mv ./*.fastq.gz ./unTrimmed_reads

  echo 'Total number of samples: '${#list[@]}

  assignThreads
  myNum=0

  for tmpFile in tmpList*tmp
  do

    # From the tmpList_XXX.tmp file (which is a chunk of the original list) 
    # put every line as an element of a bash array.
    myList=(`cat ${tmpFile} | tr '\n' ' '`)

    echo -e 'File '${tmpFile}' contains the following samples:\n'${myList[@]}
    myNum=`expr ${myNum} + ${#myList[@]}`; echo -e 'No. of samples assigned: '${myNum}'\n'

   ( for p in ${myList[@]}
    do sleep 5

      if [[ ${SEorPE} == 3 ]]; # For single-end sequencing
      then

        echo $(date) 'Trimming reads from '${p}
        java -jar ${Trimmomatic} SE -threads 1 -quiet \
        ${WD}/reads/unTrimmed_reads/${p}.fastq.gz ${WD}/reads/${p}.fastq.gz \ 
        ILLUMINACLIP:${adapters}:2:20:9:2 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:5 MINLEN:36

      else; # For paired-end sequencing

        echo $(date) 'Trimming reads from '${p}
        java -jar ${Trimmomatic} PE -threads 1 -quiet \
        ${WD}/reads/unTrimmed_reads/${p}_1.fastq.gz ${WD}/reads/unTrimmed_reads/${p}_2.fastq.gz \
        ${WD}/reads/${p}_1.fastq.gz ${WD}/reads/${p}_U1.fastq.gz \
        ${WD}/reads/${p}_2.fastq.gz ${WD}/reads/${p}_U2.fastq.gz \
        ILLUMINACLIP:${adapters}:2:20:9:2 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:5 MINLEN:36

        # Concatenate unpaired reads in a single file for every sample
        cat ${WD}/reads/${p}_U1.fastq.gz ${WD}/reads/${p}_U2.fastq.gz > ${WD}/reads/${p}_U.fastq.gz
        rm ${WD}/reads/${p}_U1.fastq.gz ${WD}/reads/${p}_U2.fastq.gz

      fi

    done ) &

  done
  wait

  rm *tmp

  ################
  doWeHaveReads
  ################

  # Remove untrimmed reads. Comment this line to avoid this behavior.
  rm -rf unTrimmed_reads

  echo -e '\nTrimming on  '${runName}' files seems to be completed'$(date)'\n'

else

  echo -e '\nTrimming will not be executed this time\n'$(date)'\n'

fi


#### ------------------------------ ####
#### ----------- Mapping ---------- ####


if [[ ${TASKS} == *M* ]]
then

  if [ ! -d ${WD}/mapping ]; then mkdir ${WD}/mapping ; echo 'The mapping path '${WD}'/mapping was created'; fi
  cd ${WD}/mapping

  echo -e '\nStarting mapping on '${runName}' files\n' $(date)'\n'

  assignThreads
  myNum=0

  for tmpFile in tmpList*tmp
  do

    # From the tmpList_XXX.tmp file (which is a chunk of the original list) 
    # put every line as an element of a bash array.
    myList=(`cat ${tmpFile} | tr '\n' ' '`)

    echo -e 'File '${tmpFile}' contains the following samples:\n'${myList[@]}
    myNum=`expr ${myNum} + ${#myList[@]}`; echo -e 'No. of samples assigned: '${myNum}'\n'

  ( for p in ${myList[@]}
    do sleep 5

    # Map the reads and sort the alignment

      echo $(date) 'Mapping reads from '${p}
      mkdir ${p}_tmpdir

      if [[ ${SEorPE} == 3 ]]; # For single-end sequencing
      then

        ${BOWTIE2} --rg-id ${p} --rg SM:${p} --rg PL:ILLUMINA -t -x \
        ${REF} -U ${WD}/reads/${p}.fastq.gz  2> ${p}_bowtie2.log | \
        java -Xmx3g -jar ${PICARD} SortSam MAX_RECORDS_IN_RAM=1000000 SO=coordinate CREATE_INDEX=true\
        TMP_DIR=${p}_tmpdir I=/dev/stdin O=${p}_bowtie2_sorted.bam >& ${p}_bowtie2_sort.log

      else; # For paired-end sequencing

        # Check that the unpaired file of reads for every sample (which was produced after trimming) exists.
        # If so, put the bowtie2 -U flag with its full location into the 'unpairandtrimmed' variable.
        # Notice this variable is called in the following line. If the file does not exist, 
        # the variable keeps empty and then the -U argument is not called in bowtie2.
        test -e ${WD}/reads/${p}_U.fastq.gz && unpairandtrimmed='-U '${WD}'/reads/'${p}'_U.fastq.gz' || unpairandtrimmed=''

        ${BOWTIE2} --rg-id ${p} --rg SM:${p} --rg PL:ILLUMINA -I 0 -X 1000 -t -x \
        ${REF} -1 ${WD}/reads/${p}_1.fastq.gz -2 ${WD}/reads/${p}_2.fastq.gz ${unpairandtrimmed}  2> ${p}_bowtie2.log | \
        java -Xmx3g -jar ${PICARD} SortSam MAX_RECORDS_IN_RAM=1000000 SO=coordinate CREATE_INDEX=true \
        TMP_DIR=${p}_tmpdir I=/dev/stdin O=${p}_bowtie2_sorted.bam >& ${p}_bowtie2_sort.log

      fi

      rm -rf ${p}_tmpdir

    # Calculate statistics from the alignments file

      echo $(date) 'Calculating statistics for '${p}

      java -Xmx3g -jar ${NGSEP} QualStats ${REF} ${p}_bowtie2_sorted.bam >& ${p}_bowtie2_readpos.stats
      java -Xmx3g -jar ${NGSEP} CoverageStats ${p}_bowtie2_sorted.bam ${p}_bowtie2_coverage.stats >& ${p}_bowtie2_coverage.log

      echo $(date) ${p}' is DONE'

    done ) &

  done
  wait

  rm *tmp

  ################
  # Check mapping failures
  numErrors=0
  for index in ${!list[@]}
  do if [[ ! `tail -1 ${list[${index}]}_bowtie2_coverage.stats` == *More* ]]
    then numErrors=`expr ${numErrors} + 1`
    echo 'Error: Mapping for sample '${list[${index}]}' failed at some point !!'
  fi; done
  if [[ ${numErrors} > 0 ]]; then echo 'Error: Mapping failed for '${numErrors}' samples'; exit 1; fi
  ################

  echo -e '\nMapping on '${runName}' files seems to be completed\n'$(date)'\n'

else

  echo -e '\nMapping will not be executed this time\n'$(date)'\n'

fi


#### ------------------------------ ####
#### ----- Variant Discovery ------ ####


if [[ ${TASKS} == *V* ]]
then

  if [ ! -d ${WD}/mapping ]; then echo 'The mapping path '${WD}'/mapping does not exist !'; exit 1; fi
  cd ${WD}/mapping

  echo -e '\nStarting variant discovery on '${runName}' files\n'$(date)'\n'

  assignThreads
  myNum=0

  for tmpFile in tmpList*tmp
  do

    # From the tmpList_XXX.tmp file (which is a chunk of the original list) 
    # put every line as an element of a bash array.
    myList=(`cat ${tmpFile} | tr '\n' ' '`)

    echo -e 'File '${tmpFile}' contains the following samples:\n'${myList[@]}
    myNum=`expr ${myNum} + ${#myList[@]}`; echo -e 'No. of samples assigned: '${myNum}'\n'

  ( for p in ${myList[@]}
    do sleep 5

      echo $(date) 'Variant discovery in '${p}

      java -Xmx3g -jar ${NGSEP} FindVariants -h 0.0001 -maxBaseQS 30 -minQuality 40 \
      -noRep -noRD -noRP -maxAlnsPerStartPos 100 -ignore5 ${i5} -ignore3 ${i3} -sampleId ${p} \
      -knownSTRs ${STRs} ${REF} ${p}_bowtie2_sorted.bam ${p}_bowtie2_NGSEP >& ${p}_bowtie2_NGSEP.log
      bgzip ${p}_bowtie2_NGSEP.vcf

    done ) &

  done
  wait

  rm *tmp

  ################
  # Check variant discovery failures
  numErrors=0
  for index in ${!list[@]}
  do if [[ ! `tail -1 ${list[${index}]}_bowtie2_NGSEP.log` == *Completed* ]]
    then numErrors=`expr ${numErrors} + 1`
    echo 'Error: Variant discovery for sample '${list[${index}]}' failed at some point !!'
  fi; done
  if [[ ${numErrors} > 0 ]]; then echo 'Error: Variant discovery failed for '${numErrors}' samples'; exit 1; fi
  ################

  echo -e '\nVariant-discovery on '${runName}' files seems to be completed\n'$(date)'\n'

else

  echo -e '\nThe FindVariants process will not be executed this time\n'$(date)'\n'

fi
