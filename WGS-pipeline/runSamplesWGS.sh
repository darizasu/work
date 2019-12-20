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
WD=/bioinfo1/projects/bean/WGS/v2.1

# Give a name for this run
runName=WGSnew

# The number of subprocesses you want to run. It depends on the number of available cores.
numThreads=16

# Specify the task(s) you want to perform. Include only the initial capital letter in a single string.
# It can include 'I'nsert-Length test, 'T'rimming, 'M'apping, 'V'ariant-Discovery.
# To run all these analyses: ./runSamples_WGS.sh 'ITMV'
TASKS=$1

# This file must (MUST) be located at ${WD}/reads. This is a fasta file containing adapter 
# sequences to be removed from the deconvoluted reads. Check Trimmomatic manual for more info.
# In case you do not have adapter sequences to be removed, leave this parameter empty.
adapters=/home/dariza/adaptersGBS.fa

# This file must (MUST) contain its full path, check description above.
samples2test=/bioinfo1/projects/bean/WGS/v2.1/samples_WGS.txt

  # Path to Software used

NGSEP=/home/dariza/software/NGSEP/NGSEPcore_3.3.0.jar
BOWTIE2=/home/dariza/bin/bowtie2
PICARD=/home/dariza/software/picard/picard.jar
Trimmomatic=/bioinfo1/software/Trimmomatic-0.36/trimmomatic-0.36.jar
Rscript=/home/dariza/bin/Rscript
BGZIP=/usr/bin/bgzip

  # Reference genome files

REF=/bioinfo1/references/bean/G19833/v2.1/bowtie2/Pvulgaris_442_v2.0.fa
STRs=/bioinfo1/references/bean/G19833/v2.1/STRs/Pvulgaris_v2_strs.list


#### ---------------------------------------------------------------------- ####
#### ----------- DO NOT MODIFY ANYTHING FROM THIS POINT FORWARD ----------- ####
#### ---------------------------------------------------------------------- ####


# Make sure you DON'T have ASCII text with CRLF line terminators
dos2unix ${samples2test}

echo -e '\nThis run was executed by:  '$(whoami)'\n'

list=(`grep -v '#' ${samples2test} | cut -f 1 | tr '\n' ' '`)
echo -e 'This run contains the following samples:\n'${list[@]}'\n'

# The whole list of samples is divided in 'nThreads' sublists.
# Every sublist is run in the background. It is an individual file 
# called tmpList_XXX.tmp that contains chunks of the original list of samples.
# It won't continue until all the sublists have finished.
function assignThreads {

  # 1st argument is a bash array with sample names
  arr=("$@")
  # 2nd argument is the number of threads to be used in the machine
  nt=${arr[-1]}
  unset 'arr[${#arr[@]}-1]'

  for i in ${!arr[@]}
  do

    samplesPerList=`expr ${i} % ${nt}`
    i=`expr ${i} + 1`
    grep -v '#' ${samples2test} | sed -n -e ${i}p >> tmpList_${samplesPerList}.tmp
  
  done
}

# Verify successful operation using a last keyword in log files
function checkLastWord {

  nErr=0
  samN=()

  # 1st argument is a bash array with sample names
  arr=("$@")
  # 2nd argument is the last word to be searched in a file
  word=${arr[-1]}
  unset 'arr[${#arr[@]}-1]'
  # 3rd argument is the extension of the file to be searched
  ext=${arr[-1]}
  unset 'arr[${#arr[@]}-1]'

  for s in ${arr[@]}
  do
    if [[ ! `tail -1 ${s}${ext}` == *${word}* ]]
    then
        nErr=`expr ${nErr} + 1`
        samN=(${samN[@]} ${s})
    fi
  done

  test ${nErr} -gt 0 && \
  echo 'This task failed for '${nErr}' samples:' && \
  echo ${samN[@]} && \
  exit 1

}

# Print the distribution of samples per thread
function printThreads {

  tmp=$1
  lst=(`cat ${tmp} | cut -f1 | tr '\n' ' '`)

    echo -e 'File '${tmp}' contains the following samples:\n'${lst[@]}
    myNum=`expr ${myNum} + ${#lst[@]}`
    echo -e 'No. of samples assigned: '${myNum}'\n'

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

  assignThreads ${list[@]} ${numThreads}
  myNum=0
  sed -i -e 's/^/#/' ${samples2test}

  for tmpFile in tmpList*tmp
  do

    # From the tmpList_XXX.tmp file (which is a chunk of the original list) 
    # put every line as an element of a bash array.
    myList=(`grep -v '#' ${tmpFile} | cut -f 1 | tr '\n' ' '`)

    test ! ${printIt} && printThreads ${tmpFile}

   ( for p in ${myList[@]}
    do sleep 1

      zcat ${WD}/reads/${p}_1.fastq.gz | head -n 4000000 > ${p}_testInsert_1.fastq
      zcat ${WD}/reads/${p}_2.fastq.gz | head -n 4000000 > ${p}_testInsert_2.fastq

      mkdir ${p}_tmpdir

      echo $(date) 'Aligning reads from '${p}

      ${BOWTIE2} --rg-id ${p} --rg SM:${p} --rg PL:ILLUMINA -I 0 -X 1200 -t -x ${REF} \
      -1 ${p}_testInsert_1.fastq -2 ${p}_testInsert_2.fastq 2> ${p}_testInsert.log | \
      java -Xmx3g -jar ${PICARD} SortSam MAX_RECORDS_IN_RAM=1000000 SO=coordinate \
      CREATE_INDEX=true TMP_DIR=${p}_tmpdir I=/dev/stdin O=${p}_testInsert_sorted.bam \
      > ${p}_InsertSizeMetrics.log 2>&1

      echo $(date) 'Collecting insert size metrics for '${p}

      java -jar ${PICARD} CollectInsertSizeMetrics HISTOGRAM_FILE=${p}_InsertSizeMetrics.pdf \
      INPUT=${p}_testInsert_sorted.bam OUTPUT=${p}_InsertSizeMetrics.txt \
      >> ${p}_InsertSizeMetrics.log 2>&1

      # Calculate Max Insert size using +/- 4 MADs from the median

      grep -A 2 METRICS ${p}_InsertSizeMetrics.txt | tail -n 2 | \
        python -c "import sys; print('\n'.join('\t'.join(c) for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip()))))" < \
        /dev/stdin | grep -e MEDIAN | cut -f2 > ${p}.metrics

        med=`head -n 1 ${p}.metrics`
        mad=`tail -n 1 ${p}.metrics`

        maxIns=`expr ${med} + 4 \* ${mad}`
        minIns=`expr ${med} - 4 \* ${mad}`

        test ${minIns} -lt 0 && minIns=0

        echo -e ${p}'\t'${minIns}'\t'${maxIns} >> ${samples2test}

     # remove temporal files

      rm -f ${p}_testInsert_*.fastq ${p}_testInsert_sorted.ba* ${p}_testInsert.log
      rm -r ${p}_tmpdir ${p}.metrics

    done ) &

  done
  wait

  ################
  checkLastWord ${list[@]} _InsertSizeMetrics.log totalMemory
  ################

  rm tmpList*tmp *InsertSizeMetrics.log

  printIt="stop"

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

  list=(`grep -v '#' ${samples2test} | cut -f 1 | tr '\n' ' '`)

  for samp in ${list[@]}; do mv ./${samp}*.fastq.gz ./unTrimmed_reads; done

  echo 'Total number of samples: '${#list[@]}

  assignThreads ${list[@]} ${numThreads}
  myNum=0

  for tmpFile in tmpList*tmp
  do

    # From the tmpList_XXX.tmp file (which is a chunk of the original list) 
    # put every line as an element of a bash array.
    myList=(`cat ${tmpFile} | cut -f 1 | tr '\n' ' '`)

    test ! ${printIt} && printThreads ${tmpFile}

  ( for p in ${myList[@]}
    do

      sleep 1

      echo $(date) 'Trimming reads from '${p}
      java -jar ${Trimmomatic} PE -threads 1 ${WD}/reads/unTrimmed_reads/${p}_1.fastq.gz \
      ${WD}/reads/unTrimmed_reads/${p}_2.fastq.gz ${WD}/reads/${p}_1.fastq.gz ${WD}/reads/${p}_U1.fastq.gz \
      ${WD}/reads/${p}_2.fastq.gz ${WD}/reads/${p}_U2.fastq.gz \
      ILLUMINACLIP:${adapters}:2:20:9:2:true LEADING:5 TRAILING:5 SLIDINGWINDOW:4:5 MINLEN:36 > ${p}_trimmomatic.log 2>&1

      cat ${WD}/reads/${p}_U1.fastq.gz ${WD}/reads/${p}_U2.fastq.gz > ${WD}/reads/${p}_U.fastq.gz

    done ) &

  done
  wait

  ################
  checkLastWord ${list[@]} _trimmomatic.log Completed
  ################

  rm *tmp #*trimmomatic.log

  printIt="stop"

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

  list=(`grep -v '#' ${samples2test} | cut -f 1 | tr '\n' ' '`)
  echo 'Total number of samples: '${#list[@]}

  assignThreads ${list[@]} ${numThreads}
  myNum=0
  sed -i -e 's/^/#/' ${samples2test}

  for tmpFile in tmpList*tmp
  do

    # From the tmpList_XXX.tmp file (which is a chunk of the original list) 
    # put every column as an element of a bash array.
    myList=(`cut -f 1 ${tmpFile} | tr '\n' ' '`)
    I=(`cut -f 2 ${tmpFile} | tr '\n' ' '`)
    X=(`cut -f 3 ${tmpFile} | tr '\n' ' '`)

    test ! ${printIt} && printThreads ${tmpFile}

   ( for index in ${!myList[@]}
    do sleep 1

   # Map the reads and sort the alignment using BOWTIE2

      echo $(date) 'Mapping reads from '${myList[${index}]}
      mkdir ${myList[${index}]}_tmpdir

      test -e ${WD}/reads/${myList[${index}]}_U.fastq.gz \
      && unpairandtrimmed='-U '${WD}'/reads/'${myList[${index}]}'_U.fastq.gz' || unpairandtrimmed=''

      ${BOWTIE2} --rg-id ${myList[${index}]} --rg SM:${myList[${index}]} --rg PL:ILLUMINA -k 3 \
      -I ${I[${index}]} -X ${X[${index}]} -t -x ${REF} -1 ${WD}/reads/${myList[${index}]}_1.fastq.gz \
      -2 ${WD}/reads/${myList[${index}]}_2.fastq.gz -p 1 ${unpairandtrimmed}  2> ${myList[${index}]}_bowtie2.log | \
      java -Xmx30g -jar ${PICARD} SortSam MAX_RECORDS_IN_RAM=1000000 SO=coordinate CREATE_INDEX=true \
      TMP_DIR=${myList[${index}]}_tmpdir I=/dev/stdin O=${myList[${index}]}_bowtie2_sorted.bam \
      >& ${myList[${index}]}_bowtie2_sort.log

      rm -rf ${myList[${index}]}_tmpdir
      rm ${WD}/reads/${myList[${index}]}_1.fastq.gz \
        ${WD}/reads/${myList[${index}]}_2.fastq.gz \
        ${WD}/reads/${myList[${index}]}_U.fastq.gz

    # Calculate statistics from the alignments file

      echo $(date) 'Calculating statistics for '${myList[${index}]}

      java -Xmx5g -jar ${NGSEP} QualStats ${REF} ${myList[${index}]}_bowtie2_sorted.bam \
      >& ${myList[${index}]}_bowtie2_readpos.stats
      java -Xmx5g -jar ${NGSEP} CoverageStats ${myList[${index}]}_bowtie2_sorted.bam \
      ${myList[${index}]}_bowtie2_coverage.stats >& ${myList[${index}]}_bowtie2_coverage.log

      ${Rscript} -e "source('https://raw.githubusercontent.com/darizasu/work/master/WGS-pipeline/readPosThreshold.R')" \
        ${myList[${index}]} ${I[${index}]} ${X[${index}]} >> ${samples2test}

      echo $(date) ${myList[${index}]}' is DONE'

    done ) &

  done
  wait

  rm *tmp

  printIt="stop"

  ################
  checkLastWord ${list[@]} _bowtie2_coverage.stats More
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

  list=(`grep -v '#' ${samples2test} | cut -f 1 | tr '\n' ' '`)
  echo 'Total number of samples: '${#list[@]}

  assignThreads ${list[@]} ${numThreads}
  myNum=0

  for tmpFile in tmpList*tmp
  do

    # From the tmpList_XXX.tmp file (which is a chunk of the original list) 
    # put every column as an element of a bash array.
    myList=(`cut -f 1 ${tmpFile} | tr '\n' ' '`)
    Ifive=(`cut -f 4 ${tmpFile} | tr '\n' ' '`)
    Ithree=(`cut -f 5 ${tmpFile} | tr '\n' ' '`)

    test ! ${printIt} && printThreads ${tmpFile}

   ( for index in ${!myList[@]}
    do sleep 1

      echo $(date) 'Variant discovery in '${myList[${index}]}

      java -Xmx15g -jar ${NGSEP} FindVariants -h 0.0001 -noRep -maxBaseQS 30 -minQuality 40 \
      -maxAlnsPerStartPos 2 -ignore5 ${Ifive[${index}]} -ignore3 ${Ithree[${index}]} \
      -sampleId ${myList[${index}]} -knownSTRs ${STRs} ${REF} ${myList[${index}]}_bowtie2_sorted.bam \
      ${myList[${index}]}_bowtie2_NGSEP >& ${myList[${index}]}_bowtie2_NGSEP.log
      ${BGZIP} ${myList[${index}]}_bowtie2_NGSEP.vcf

    done ) &

  done
  wait

  rm *tmp

  printIt="stop"

  ################
  checkLastWord ${list[@]} _bowtie2_NGSEP.log Completed
  ################

  echo -e '\nVariant-discovery on '${runName}' files seems to be completed\n'$(date)'\n'

else

  echo -e '\nThe FindVariants process will not be executed this time\n'$(date)'\n'

fi
