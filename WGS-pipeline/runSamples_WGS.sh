#!/bin/bash

# This wrapper script is intended to run Mapping (using Bowtie2) and FindVariants
# (using NGSEP) for a list of WGS samples.
# Please follow these instructions to run this script properly.

# The parameter 'samples2test' specifies the full path and filename to
# the file that contains parameters for Mapping and FindVariants.
# From the begining, this file is created, used, then modified the 1st time,
# used, then modified the 2nd time, and used for the last time.
# In the 1st step (see below), 'samples2test' contains one line for each sample_id
# In the 2nd step (see below), 'samples2test' contains three tab-separated columns:
# sample_id[tab]minins[tab]maxins.
# In the 3rd step (see below), 'samples2test' contains five tab-separated columns:
# sample_id[tab]minins[tab]maxins[tab]I5prime[tab]I3prime
# Please follow these steps.
#
  # 1) Create a txt file containing a sample_id per line.
  # Comment lines with '#' are accepted.
  # Specify the name of this file in the 'samples2test' parameter as described above.
  # Execute the 1st run of 'runSamples_WGS.sh' specifying the task 'I'
  # (just for InsertLength test)
  # It takes the first 1.000.000 reads for each sample's reads and maps them to
  # the reference genome.
  # Then it runs Picard's CollectInsertSizeMetrics module to get the insert size
  # characteristics for each sample's sequencing run
  # After running this module, a pdf histogram and a txt files are produced showing
  # the insert size distribution.
  # You have to decide the minimum and maximum insert size length for each sample by 
  # checking the points where the slope starts to increase and decrease rapidly
  # (minimum and maximum insert size length respectively).
  # minins is the minimum fragment length for a valid paired-end alignment.
  # maxins is the maximum fragment length for a valid paired-end alignment.
  # Check Bowtie2 manual page for more info.
  # Once you have decided the minins and maxins parameters, specify them in
  # 'samples2test' file for each sample_id as the second and third tab-separated columns respectively.
#
  # 2) Execute the 2nd run of 'runSamples_WGS.sh' specifying the tasks 'TM'
  # (just for trimming and mapping)
  # Once correctly completed, the files *_bowtie2_readpos.stats are produced by
  # using NGSEP-QualStats.
  # Plot these results. You have to decide the I5prime and I3prime parameters by
  # looking at how many consecutive positions in the 5' end and the 3' end have a
  # high sequencing error bias (guide yourself by the behavior of the slope).
  # I5prime is to ignore this many base pairs from the 5' end of the reads.
  # I3prime is to ignore this many base pairs from the 3' end of the reads,
  # check NGSEP manual page. Once you have decided the I5prime and I3prime parameters,
  # specify them in 'samples2test' for each sample_id as the fourth and fifth tab-separated columns respectively.
#
  # 3) Exectute the 3rd run of 'runMapping_WGS.sh' specifying the task 'V'

  # Key variables to specify

# This is the working directory full path. It must contain a subdirectory 'reads'
# where the raw reads for the samples to be processed are located.
WD=/bioinfo2/projects/bean_tmp/WGS_to_v2.1

# Give a name for this run
runName=WGS

# The number of subprocesses you want to run. It depends on the number of available cores.
numThreads=5

# Specify the task(s) you want to perform. Include only the initial capital letter in a single string.
# It can include 'I'nsert-Length test, 'T'rimming, 'M'apping, 'V'ariant-Discovery.
# To run all these analyses: ./runSamples_WGS.sh 'ITMV'
TASKS=$1

# This file must (MUST) be located at ${WD}/reads. This is a fasta file containing adapter 
# sequences to be removed from the deconvoluted reads. Check Trimmomatic manual for more info.
# In case you do not have adapter sequences to be removed, leave this parameter empty.
adapters=

# This file must (MUST) contain its full path, check description above.
samples2test=/bioinfo2/projects/bean_tmp/WGS_to_v2.1/mapping/samples2test.txt; 

  # Path to Software used

NGSEP=/data/software/NGSEPcore_3.0.2.jar
BOWTIE2=/data/software/bowtie2-2.2.9/bowtie2
PICARD=/data/software/picard-tools-1.140/picard.jar
Trimmomatic=/home/dariza/software/Trimmomatic-0.36

  # Reference genome files

REF=/data/references/bean/v2.1/bowtie2/Pvulgaris_442_v2.0.fa
STRs=/data/references/bean/v2.1/strs/Pvulgaris_v2_strs.list


#### ---------------------------------------------------------------------- ####
#### ----------- DO NOT MODIFY ANYTHING FROM THIS POINT FORWARD ----------- ####
#### ---------------------------------------------------------------------- ####


echo -e '\nThis run was executed by:  '$(whoami)'\n'

list=(`grep -v '#' ${samples2test} | cut -f 1 | tr '\n' ' '`)
echo -e 'This run contains the following samples:\n'${list[@]}'\n'

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


#### ------------------------------ ####
#### -------- Insert Length ------- ####


if [[ ${TASKS} == *I* ]]
then

  if [ ! -d ${WD}/reads ]; then echo 'The directory '${WD}'/reads does not exist !'; exit 1; fi

  if [ ! -d ${WD}/mapping ]
    then mkdir ${WD}/mapping ; echo 'The mapping path '${WD}'/mapping was created'
  fi
  
  cd ${WD}/mapping

  echo -e '\nStarting Insert Length test on '${runName}' files\n'$(date)'\n'
  echo 'Total number of samples: '${#list[@]}

  assignThreads
  myNum=0

  for tmpFile in tmpList*tmp
  do

    # From the tmpList_XXX.tmp file (which is a chunk of the original list) 
    # put every line as an element of a bash array.
    myList=(`grep -v '#' ${tmpFile} | cut -f 1 | tr '\n' ' '`)
    echo -e 'File '${tmpFile}' contains the following samples:\n'${myList[@]}
    myNum=`expr ${myNum} + ${#myList[@]}`; echo -e 'No. of samples assigned: '${myNum}'\n'

   ( for p in ${myList[@]}
    do sleep 5

      zcat ${WD}/reads/${p}_1.fastq.gz | head -n 4000000 > ${p}_testInsert_1.fastq;
      zcat ${WD}/reads/${p}_2.fastq.gz | head -n 4000000 > ${p}_testInsert_2.fastq;

      mkdir ${p}_tmpdir

      echo $(date) 'Aligning reads from '${p}

      ${BOWTIE2} --rg-id ${p} --rg SM:${p} --rg PL:ILLUMINA -I 0 -X 1000 -t -x ${REF} \
      -1 ${p}_testInsert_1.fastq -2 ${p}_testInsert_2.fastq 2> ${p}_testInsert.log | \
      java -Xmx3g -jar ${PICARD} SortSam MAX_RECORDS_IN_RAM=1000000 SO=coordinate \
      CREATE_INDEX=true TMP_DIR=${p}_tmpdir I=/dev/stdin O=${p}_testInsert_sorted.bam \
      > ${p}_InsertSizeMetrics.log 2>&1

      echo $(date) 'Collecting insert size metrics for '${p}

      java -jar ${PICARD} CollectInsertSizeMetrics HISTOGRAM_FILE=${p}_InsertSizeMetrics.pdf \
      INPUT=${p}_testInsert_sorted.bam OUTPUT=${p}_InsertSizeMetrics.txt \
      >> ${p}_InsertSizeMetrics.log 2>&1

     # remove temporal files

      rm -f ${p}_testInsert_*.fastq ${p}_testInsert_sorted.ba* ${p}_testInsert.log
      rm -r ${p}_tmpdir

    done ) &

  done
  wait

  ################
  for index in ${!list[@]}; do
    if [[ ! `tail -1 ${list[${index}]}_InsertSizeMetrics.log` == *totalMemory* ]]
      then echo "Error: Test insert lenght failed for "${list[${index}]}" at some point !!"; exit 1
    fi
  done
  ################

  rm tmpList*tmp ${p}_InsertSizeMetrics.log

  echo -e '\nInsert Length test on '${runName}' files seems to be completed\n'$(date)'\n'

else

  echo -e '\nTest insert lenght will not be executed this time\n'$(date)'\n'

fi


#### ------------------------------ ####
#### ----------- Trimming --------- ####


if [[ ${TASKS} == *T* ]]
then

  if [ ! -d ${WD}/reads ]; then echo 'The reads path '${WD}'/reads does not exist !'; exit 1; fi
  cd ${WD}/reads

  echo -e '\nStarting trimming on '${runName}' files'$(date)'\n'
  mkdir unTrimmed_reads

  for samp in ${list[@]}; do mv ./${samp}*.fastq.gz ./unTrimmed_reads; done

  echo 'Total number of samples: '${#list[@]}

  assignThreads
  myNum=0

  for tmpFile in tmpList*tmp
  do

    # From the tmpList_XXX.tmp file (which is a chunk of the original list) 
    # put every line as an element of a bash array.
    myList=(`cat ${tmpFile} | cut -f 1 | tr '\n' ' '`)

    echo 'File '${tmpFile}' contains the following samples: '; echo ${myList[@]}
    myNum=`expr ${myNum} + ${#myList[@]}`; echo -e 'No. of samples assigned: '${myNum}'\n'

  ( for p in ${myList[@]}
    do

      sleep 5

      echo $(date) 'Trimming reads from '${p}
      java -jar ${Trimmomatic}/trimmomatic-0.36.jar PE -threads 1 ${WD}/reads/unTrimmed_reads/${p}_1.fastq.gz \
      ${WD}/reads/unTrimmed_reads/${p}_2.fastq.gz ${WD}/reads/${p}_1.fastq.gz ${WD}/reads/${p}_U1.fastq.gz \
      ${WD}/reads/${p}_2.fastq.gz ${WD}/reads/${p}_U2.fastq.gz \
      ILLUMINACLIP:${adapters}:2:20:9:2 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:5 MINLEN:36 > ${p}_trimmomatic.log 2>&1

      cat ${WD}/reads/${p}_U1.fastq.gz ${WD}/reads/${p}_U2.fastq.gz > ${WD}/reads/${p}_U.fastq.gz

    done ) &

  done
  wait

  ################
  for index in ${!list[@]}; do
    if [[ ! `tail -1 ${list[${index}]}_trimmomatic.log` == *Completed* ]]
      then echo "Error: There are no trimmed reads for "${list[${index}]}; exit 1
    fi
  done
  ################

  rm *tmp *trimmomatic.log

  # Remove untrimmed reads. Comment this line to avoid this behavior.
#  rm -rf unTrimmed_reads

  echo -e '\nTrimming on  '${runName}' files seems to be completed\n'$(date)'\n'

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

  echo 'Total number of samples: '${#list[@]}

  assignThreads
  myNum=0

  for tmpFile in tmpList*tmp
  do

    # From the tmpList_XXX.tmp file (which is a chunk of the original list) 
    # put every column as an element of a bash array.
    myList=(`cut -f 1 ${tmpFile} | tr '\n' ' '`)
    I=(`cut -f 2 ${tmpFile} | tr '\n' ' '`)
    X=(`cut -f 3 ${tmpFile} | tr '\n' ' '`)

    echo -e 'File '${tmpFile}' contains the following samples:\n'${myList[@]}
    myNum=`expr ${myNum} + ${#myList[@]}`; echo -e 'No. of samples assigned: '${myNum}'\n'

   ( for index in ${!myList[@]}
    do sleep 5

   # Map the reads and sort the alignment using BOWTIE2

      echo $(date) 'Mapping reads from '${myList[${index}]}
      mkdir ${myList[${index}]}_tmpdir

      test -e ${WD}/reads/${myList[${index}]}_U.fastq.gz \
      && unpairandtrimmed='-U '${WD}'/reads/'${myList[${index}]}'_U.fastq.gz' || unpairandtrimmed=''

      ${BOWTIE2} --rg-id ${myList[${index}]} --rg SM:${myList[${index}]} --rg PL:ILLUMINA -k 3 \
      -I ${I[${index}]} -X ${X[${index}]} -t -x ${REF} -1 ${WD}/reads/${myList[${index}]}_1.fastq.gz \
      -2 ${WD}/reads/${myList[${index}]}_2.fastq.gz ${unpairandtrimmed}  2> ${myList[${index}]}_bowtie2.log | \
      java -Xmx10g -jar ${PICARD} SortSam MAX_RECORDS_IN_RAM=1000000 SO=coordinate CREATE_INDEX=true \
      TMP_DIR=${myList[${index}]}_tmpdir I=/dev/stdin O=${myList[${index}]}_bowtie2_sorted.bam \
      >& ${myList[${index}]}_bowtie2_sort.log

      rm -rf ${myList[${index}]}_tmpdir

    # Calculate statistics from the alignments file

      echo $(date) 'Calculating statistics for '${myList[${index}]}

      java -Xmx5g -jar ${NGSEP} QualStats ${REF} ${myList[${index}]}_bowtie2_sorted.bam \
      >& ${myList[${index}]}_bowtie2_readpos.stats
      java -Xmx5g -jar ${NGSEP} CoverageStats ${myList[${index}]}_bowtie2_sorted.bam \
      ${myList[${index}]}_bowtie2_coverage.stats >& ${myList[${index}]}_bowtie2_coverage.log

      echo $(date) ${myList[${index}]}' is DONE'

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

  echo 'Total number of samples: '${#list[@]}

  assignThreads
  myNum=0

  for tmpFile in tmpList*tmp
  do

    # From the tmpList_XXX.tmp file (which is a chunk of the original list) 
    # put every column as an element of a bash array.
    myList=(`cut -f 1 ${tmpFile} | tr '\n' ' '`)
    Ifive=(`cut -f 4 ${tmpFile} | tr '\n' ' '`)
    Ithree=(`cut -f 5 ${tmpFile} | tr '\n' ' '`)

    echo 'File '${tmpFile}' contains the following samples: '; echo ${myList[@]}
    myNum=`expr ${myNum} + ${#myList[@]}`; echo -e 'No. of samples assigned: '${myNum}'\n'

   ( for index in ${!myList[@]}
    do sleep 5

      echo $(date) 'Variant discovery in '${myList[${index}]}

      java -Xmx15g -jar ${NGSEP} FindVariants -h 0.0001 -maxBaseQS 30 -minQuality 40 \
      -maxAlnsPerStartPos 2 -ignore5 ${Ifive[${index}]} -ignore3 ${Ithree[${index}]} \
      -sampleId ${myList[${index}]} -knownSTRs ${STRs} ${REF} ${myList[${index}]}_bowtie2_sorted.bam \
      ${myList[${index}]}_bowtie2_NGSEP >& ${myList[${index}]}_bowtie2_NGSEP.log
      bgzip ${myList[${index}]}_bowtie2_NGSEP.vcf

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
