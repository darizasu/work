#!/bin/bash

  # Key variables to specify

# This is the working directory full path. It should contain two directories: 
# 'reads' and 'mapping'. 'reads' must have a subdirectory called 'lane', 
# which contains the raw reads.
WD=/bioinfo1/projects/bean/GBSplates/23

# This file should be located at ${WD}/reads/lane , otherwise its path must be specified. 
# Check NGSEP Deconvolute <INDEX_FILE> parameter for more info.
INDEXFILE=/bioinfo1/projects/bean/GBSplates/23/reads/lane/barcodeMap_plate23.txt

# This file should be located at ${WD}/reads/lane , otherwise its path must be specified. 
# Check NGSEP Deconvolute -d flag for more info.
FILES2DECONV=/bioinfo1/projects/bean/GBSplates/23/reads/lane/lanes_plate23.txt

# This is your plate's name
runName=plate_23

# The number of subprocesses you want to run. It depends on the number of available cores.
numThreads=6

# Specify the task(s) you want to perform. Include only the initial capital letter in a single string.
# It can include 'D'econvolution, 'T'rimming, 'M'apping, 'V'ariant-Discovery.
# To run all these analyses: ./runPlate.sh 'DTMV'
TASKS=$1

# Ignore this many base pairs from the 5' and 3' ends in NGSEP FindVariants.
# To get this numbers, run this script with the tasks 'DTM', then run the 'calculateReadPosGBS.sh' script 
# to get sequencing error bias for the entire plate and plot those results. Then decide the i5 
# and i3 parameters and specify them in the following lines. 
# Then run again this script with the task 'V'.
i5=6
i3=12

# This file must (MUST) be located at ${WD}/reads. This is a fasta file containing adapter 
# sequences to be removed from the deconvoluted reads. Check Trimmomatic manual for more info.
adapters=adapters_21-26.fa

  # Path to Software used

NGSEP=/data/software/NGSEPcore_3.0.2.jar
BOWTIE2=/data/software/bowtie2-2.3.0/bowtie2
PICARD=/data/software/picard-tools-1.140/picard.jar
Trimmomatic=/bioinfo1/software/Trimmomatic-0.36/trimmomatic-0.36.jar

  # Reference genome files

REF=/data/references/bean/v2.1/bowtie2/Pvulgaris_442_v2.0.fa
STRs=/data/references/bean/v2.1/strs/Pvulgaris_v2_strs.list


#### ---------------------------------------------------------------------- ####
#### ----------- DO NOT MODIFY ANYTHING FROM THIS POINT FORWARD ----------- ####
#### ---------------------------------------------------------------------- ####


echo -e '\nThis run was executed by:  '$(whoami)'\n'

if [ ! -d ${WD}/reads/lane ]; then echo 'The raw reads path '${WD}'/reads/lane does not exist !'; exit 1; fi
cd ${WD}/reads/lane

# Get sequencing run type, either single-end or paired-end from 
# the 'FILES2DECONV' file based on the number of columns it contains.
SEorPE=`head -1 ${FILES2DECONV} | tr '\t' '\n' | wc -l`

test ${SEorPE} == 4 && echo -e 'This run is working with Paired-end sequencing\n' \
  || echo -e 'This run is working with Single-end sequencing\n'

# Get the list of samples from the 4th column in 'INDEXFILE' file, avoid first line. 
list=(`cut -f 4 ${INDEXFILE} | tail -n +2 | sort -u | tr '\n' ' '`)

echo -e 'This run contains the following samples:\n'${list[@]}'\n'

# The whole list of samples is divided in 'nThreads' sublists.
# Every sublist is run in the background. It is an individual file 
# called tmpList_XXX.tmp that contains chunks of the original list of samples.
# It won't continue until all the sublists have finished.
function assignThreads {
  for i in ${!list[@]}
  do  samplesPerList=`expr ${i} % ${numThreads}`
    echo ${list[${i}]} >> tmpList_${samplesPerList}.tmp
  done
}

# Check there are reads for every sample, exit otherwise.
function doWeHaveReads {
  numErrors=0
  for i in ${!list[@]}
  do if ! ls ${WD}/reads/${list[${i}]}*.fastq.gz > /dev/null 2>&1
    then numErrors=`expr ${numErrors} + 1`; echo 'Error: The sample '${list[${i}]}' does not have reads'
  fi; done
  if [[ ${numErrors} > 0 ]];
    then echo 'Error: There are '${numErrors}' samples with no reads at '${WD}'/reads/'; exit 1
  fi
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

        mv ./${p}.fastq.gz ./unTrimmed_reads

        echo $(date) 'Trimming reads from '${p}
        java -jar ${Trimmomatic} SE -threads 1 ${WD}/reads/unTrimmed_reads/${p}.fastq.gz \
        ${WD}/reads/${p}.fastq.gz ILLUMINACLIP:${adapters}:2:20:7:2 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:5 MINLEN:36 > ${p}_trimmomatic.log 2>&1

      else # For paired-end sequencing

        mv -t ./unTrimmed_reads ./${p}_1.fastq.gz ./${p}_2.fastq.gz

        echo $(date) 'Trimming reads from '${p}
        java -jar ${Trimmomatic} PE -threads 1 \
        ${WD}/reads/unTrimmed_reads/${p}_1.fastq.gz ${WD}/reads/unTrimmed_reads/${p}_2.fastq.gz \
        ${WD}/reads/${p}_1.fastq.gz ${WD}/reads/${p}_U1.fastq.gz ${WD}/reads/${p}_2.fastq.gz \
        ${WD}/reads/${p}_U2.fastq.gz ILLUMINACLIP:${adapters}:2:20:7:2 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:5 MINLEN:36 > ${p}_trimmomatic.log 2>&1

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

  for i in ${!myList[@]}; do
    if [[ ! `tail -1 ${myList[${i}]}_trimmomatic.log` == *Completed* ]]
      then echo "Error: There was a Trimmomatic error for "${myList[${i}]}; exit 1
    fi
  done

  ################

  # Remove untrimmed reads and trimmomatic logs. Comment this line to avoid this behavior.
  rm -rf unTrimmed_reads
  rm *_trimmomatic.log

  echo -e '\nTrimming on  '${runName}' files seems to be completed '$(date)'\n'

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

      else # For paired-end sequencing

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
  for i in ${!list[@]}
  do if [[ ! `tail -1 ${list[${i}]}_bowtie2_coverage.stats` == *More* ]]
    then numErrors=`expr ${numErrors} + 1`
    echo 'Error: Mapping for sample '${list[${i}]}' failed at some point !!'
  fi; done
  if [[ ${numErrors} > 0 ]]; then echo 'Error: Mapping failed for '${numErrors}' samples'; exit 1; fi
  ################

  echo -e '\n'$(date) 'Calculating Quality statistics for '${runName}

  cd ${WD}/mapping

  for col in {2..5}
  do

    touch tmpFile_${col}.tmp

    for sample in *_bowtie2_readpos.stats
    do

      head -n -3 ${sample} | cut -f${col} | paste tmpFile_${col}.tmp - > all_pl${runName}_col${col}.txt
      cp all_pl${runName}_col${col}.txt tmpFile_${col}.tmp

    done

    sed -i "s/^[ \t]*//" all_pl${runName}_col${col}.txt
    rm tmpFile_${col}.tmp

    awk '{for(i=1;i<=NF;i++) t+=$i; print t; t=0}' all_pl${runName}_col${col}.txt > sum_pl${runName}_col${col}.txt
    rm all_pl${runName}_col${col}.txt

    mv sum_pl${runName}_col${col}.txt ${WD}

  done

  cd ${WD}
  paste sum_pl${runName}_col2.txt sum_pl${runName}_col3.txt sum_pl${runName}_col4.txt \
    sum_pl${runName}_col5.txt > readPosStats_pl${runName}.tmp
  rm sum_pl*

  awk '{print $1/$3"\t"$2/$4}' readPosStats_pl${runName}.tmp > readPosStats_${runName}.txt
  rm readPosStats_pl${runName}.tmp

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
  for i in ${!list[@]}
  do if [[ ! `tail -1 ${list[${i}]}_bowtie2_NGSEP.log` == *Completed* ]]
    then numErrors=`expr ${numErrors} + 1`
    echo 'Error: Variant discovery for sample '${list[${i}]}' failed at some point !!'
  fi; done
  if [[ ${numErrors} > 0 ]]; then echo 'Error: Variant discovery failed for '${numErrors}' samples'; exit 1; fi
  ################

  echo -e '\nVariant-discovery on '${runName}' files seems to be completed\n'$(date)'\n'

else

  echo -e '\nThe FindVariants process will not be executed this time\n'$(date)'\n'

fi
