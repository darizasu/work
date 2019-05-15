#!/bin/bash

#### ---------------------------------------------------------------------- ####
#### ---------------------- KEY VARIABLES TO SPECIFY ---------------------- ####
#### ---------------------------------------------------------------------- ####

# This is the full path for the working directory .
# It must (MUST) contain the following structure inside:
   # /path/to/working/directory
   #    \_ reads
   #      \_ lane
# The subdirectory 'lane' contains the raw sequencing reads in FASTQ format
WD=/bioinfo1/projects/bean/GBS_QC/36

# Check NGSEP Deconvolute <INDEX_FILE> parameter for more information
# about the following ${INDEXFILE} file. This file should be located at ${WD}/reads/lane
# otherwise ts full path must be specified. 
INDEXFILE=/bioinfo1/projects/bean/GBS_QC/36/reads/lane/barcodeMap_plate36.txt

# Check NGSEP Deconvolute -d flag for more information about the 
# following ${FILES2DECONV} file. This file should be located at 
# ${WD}/reads/lane , otherwise its full path must be specified. 
FILES2DECONV=/bioinfo1/projects/bean/GBS_QC/36/reads/lane/lanes_plate36.txt

# This is the plate number ID.
pno=36

# The number of subprocesses you want to run.
# It depends on the number of available cores in your machine.
numThreads=10

# Specify the task(s) you want to perform. Include only the initial capital 
# letter in a single string. It can include 'D'econvolution, 'T'rimming, 
# 'M'apping, 'V'ariant-Discovery and/or
# 'Q'uality-check
# To run all these analyses: ./runPlate.sh 'DTMVQ'
TASKS=$1

# The following 'adapters' file is a FASTA file containing adapter sequences to 
# be removed from the deconvoluted reads. It is only used when you run the 
# 'T'rimming task. Check the Trimmomatic manual for more info.
adapters=

# The following parameters are used when you run the 'V'ariant-Discovery and/or
# the 'Q'uality-check tasks. They are used to ignore this many base pairs 
# from the 5' and 3' ends in NGSEP FindVariants.
# To get these numbers, run this script with the tasks 'DTM'.
# Once the mapping step is done, the file 'readPosStats_XX.txt' is created in ${WD}.
# It contains the total sequencing error bias along the read position for 
# unique-mapping reads (first column) and multi-mapping reads (second column),
# calculated for all the samples in the plate. You have to decide the i5 and i3
# parameters by looking at how many consecutive positions in the 5' and the 3' ends
# have a high sequencing error bias (guide yourself by the behavior of the slope).
# Then decide the i5 and i3 parameters and specify them in the following lines. 
# Then run again this script with the task 'V'.
i5=1
i3=10

# In case you specified the task 'Q'uality-check, 'myVariants' is a mandatory parameter.
# 'myVariants' is a comprehensive list of variants for the species of interest.
# These variants are genotyped with the mapping file that is produced in this run.
# The genotyped variants are used to compare different samples between each other.
# This is useful when identifying duplicate genotypes.
myVariants=/bioinfo1/references/bean/G19833/v2.1/variants/Variants_WGS_tpg_noGenotype.vcf.gz

  # Reference genome files

REF=/bioinfo1/references/bean/G19833/v2.1/bowtie2/Pvulgaris_442_v2.0.fa
STRs=/bioinfo1/references/bean/G19833/v2.1/STRs/Pvulgaris_v2_strs.list

  # Path to Software used

NGSEP=/home/dariza/software/NGSEP/NGSEPcore_3.3.0.jar
BOWTIE2=/data/software/bowtie2-2.3.0/bowtie2
PICARD=/data/software/picard-tools-1.140/picard.jar
Trimmomatic=/bioinfo1/software/Trimmomatic-0.36/trimmomatic-0.36.jar
BGZIP=/usr/bin/bgzip


#### ---------------------------------------------------------------------- ####
#### ----------- DO NOT MODIFY ANYTHING FROM THIS POINT FORWARD ----------- ####
#### ---------------------------------------------------------------------- ####


echo -e '\nThis run was executed by:  '$(whoami)'\n'

test ! -d ${WD}/reads/lane && \
  echo -e '\nThe raw reads path '${WD}'/reads/lane does not exist !' && \
  exit 1

cd ${WD}/reads/lane

# Make sure you DON'T have ASCII text with CRLF line terminators
dos2unix ${INDEXFILE}
dos2unix ${FILES2DECONV}

# Get sequencing run type, either single-end or paired-end from 
# the 'FILES2DECONV' file based on the number of columns it contains.
SEorPE=`head -1 ${FILES2DECONV} | tr '\t' '\n' | wc -l`

test ${SEorPE} == 4 && \
  echo -e '\nThis run is working with Paired-end sequencing\n' || \
  echo -e '\nThis run is working with Single-end sequencing\n'

# Get the list of samples from the 4th column in 'INDEXFILE' file, avoid first line. 
list=(`cut -f 4 ${INDEXFILE} | tail -n +2 | sort -u | tr '\n' ' '`)

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
    echo ${arr[${i}]} >> tmpList_${samplesPerList}_${pno}.tmp
  
  done
}

# Check there are reads for every sample, exit otherwise.
function doWeHaveReads {

  nErr=0
  samN=()

  for s in ${list[@]}
  do 

    if ! ls ${WD}/reads/${s}*.fastq.gz > /dev/null 2>&1
    then 
        nErr=`expr ${nErr} + 1`
        samN=(${samN[@]} ${s})
    fi

  done

  test ${nErr} -gt 0 && \
    echo 'Error: There are '${nErr}' samples with no reads at '${WD}'/reads/' && \
    echo ${samN[@]} && \
    samN=() &&
    exit 1
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
  lst=(`cat ${tmp} | tr '\n' ' '`)

  myNum=0
    echo -e 'File '${tmp}' contains the following samples:\n'${lst[@]}
    myNum=`expr ${myNum} + ${#lst[@]}`
    echo -e 'No. of samples assigned: '${myNum}'\n'

}


#### ------------------------------ ####
#### --------- Demultiplex -------- ####


test ! -d ${WD}/reads/lane && \
  echo 'The raw-reads path '${WD}'/reads/lane does not exist !' && \
  exit 1

cd ${WD}/reads/lane

if [[ ${TASKS} == *D* ]]
then

  echo -e '\nStarting demultiplexing on plate '${pno}' files\n'$(date)'\n'

  java -jar ${NGSEP} Demultiplex -o ${WD}/reads -d ${FILES2DECONV} ${INDEXFILE}

fi

################
doWeHaveReads
################

echo -e '\nDemultiplexing on plate '${pno}' is done\n'$(date)'\n'


#### ------------------------------ ####
#### ----------- Trimming --------- ####


if [[ ${TASKS} == *T* ]]
then

  cd ${WD}/reads

  echo -e '\nStarting trimming on plate '${pno}' '$(date)'\n'
  mkdir unTrimmed_reads

  echo 'Total number of samples: '${#list[@]}

  assignThreads ${list[@]} ${numThreads}
  
  for tmpFile in tmpList_*${pno}.tmp
  do

    # From the tmpList_XXX.tmp file (which is a chunk of the original list) 
    # put every line as an element of a bash array.
    myList=(`cat ${tmpFile} | tr '\n' ' '`)

    test ! ${printIt} && printThreads ${tmpFile}

   ( for s in ${myList[@]}
    do sleep 2

      if [[ ${SEorPE} == 3 ]] # For single-end sequencing
      then

        mv -t ./unTrimmed_reads ./${s}.fastq.gz

        echo $(date) 'Trimming reads from '${s}
        java -jar ${Trimmomatic} SE -threads 1 ${WD}/reads/unTrimmed_reads/${s}.fastq.gz \
        ${WD}/reads/${s}.fastq.gz ILLUMINACLIP:${adapters}:2:20:7:2 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:5 MINLEN:36 > ${s}_trimmomatic.log 2>&1

      else # For paired-end sequencing

        mv -t ./unTrimmed_reads ./${s}_1.fastq.gz ./${s}_2.fastq.gz

        echo $(date) 'Trimming reads from '${s}
        java -jar ${Trimmomatic} PE -threads 1 \
        ${WD}/reads/unTrimmed_reads/${s}_1.fastq.gz ${WD}/reads/unTrimmed_reads/${s}_2.fastq.gz \
        ${WD}/reads/${s}_1.fastq.gz ${WD}/reads/${s}_U1.fastq.gz ${WD}/reads/${s}_2.fastq.gz \
        ${WD}/reads/${s}_U2.fastq.gz ILLUMINACLIP:${adapters}:2:20:7:2 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:5 MINLEN:36 > ${s}_trimmomatic.log 2>&1

        # Concatenate unpaired reads in a single file for every sample
        cat ${WD}/reads/${s}_U1.fastq.gz ${WD}/reads/${s}_U2.fastq.gz > \
          ${WD}/reads/${s}_U.fastq.gz
        rm ${WD}/reads/${s}_U1.fastq.gz ${WD}/reads/${s}_U2.fastq.gz

      fi

    done ) &

  done
  wait

  rm tmpList_*${pno}.tmp

  ################
  doWeHaveReads
  checkLastWord ${list[@]} _trimmomatic.log Completed
  ################

  # Remove untrimmed reads and trimmomatic logs. Comment this line to avoid this behavior.
  rm -rf unTrimmed_reads
  rm *_trimmomatic.log

  echo -e '\nTrimming on plate '${pno}' files seems to be completed '$(date)'\n'

  printIt="stop"

else

  echo -e '\nTrimming will not be executed this time\n'$(date)'\n'

fi


#### ------------------------------ ####
#### ----------- Mapping ---------- ####


if [[ ${TASKS} == *M* ]]
then

  test ! -d ${WD}/mapping && \
    mkdir ${WD}/mapping &&
    echo 'The mapping path '${WD}'/mapping was created'

  cd ${WD}/mapping

  echo -e '\nStarting mapping on plate '${pno}' files\n' $(date)'\n'

  assignThreads ${list[@]} ${numThreads}

  for tmpFile in tmpList_*${pno}.tmp
  do

    # From the tmpList_XXX.tmp file (which is a chunk of the original list) 
    # put every line as an element of a bash array.
    myList=(`cat ${tmpFile} | tr '\n' ' '`)

    test ! ${printIt} && printThreads ${tmpFile}

  ( for s in ${myList[@]}
    do sleep 2

    # Map the reads and sort the alignment

      echo $(date) 'Mapping reads from '${s}
      mkdir ${s}_tmpdir

      if [[ ${SEorPE} == 3 ]]; # For single-end sequencing
      then

        ${BOWTIE2} --rg-id ${s} --rg SM:${s} --rg PL:ILLUMINA -t -x ${REF} \
          -U ${WD}/reads/${s}.fastq.gz  2> ${s}_bowtie2.log | \
        java -Xmx3g -jar ${PICARD} SortSam MAX_RECORDS_IN_RAM=1000000 \
          SO=coordinate CREATE_INDEX=true TMP_DIR=${s}_tmpdir I=/dev/stdin \
          O=${s}_bowtie2_sorted.bam >& ${s}_bowtie2_sort.log

      else # For paired-end sequencing

        # Check that the unpaired file of reads for every sample (which was produced after trimming) exists.
        # If so, put the bowtie2 -U flag with its full location into the 'unpairandtrimmed' variable.
        # Notice this variable is called in the following line. If the file does not exist, 
        # the variable keeps empty and then the -U argument is not called in bowtie2.
        test -e ${WD}/reads/${s}_U.fastq.gz && \
          unpairandtrimmed='-U '${WD}'/reads/'${s}'_U.fastq.gz' || \
          unpairandtrimmed=''

        ${BOWTIE2} --rg-id ${s} --rg SM:${s} --rg PL:ILLUMINA -I 0 -X 1000 -t \
          -x ${REF} -1 ${WD}/reads/${s}_1.fastq.gz -2 ${WD}/reads/${s}_2.fastq.gz \
          ${unpairandtrimmed}  2> ${s}_bowtie2.log | \
        java -Xmx3g -jar ${PICARD} SortSam MAX_RECORDS_IN_RAM=1000000 \
          SO=coordinate CREATE_INDEX=true TMP_DIR=${s}_tmpdir I=/dev/stdin \
          O=${s}_bowtie2_sorted.bam >& ${s}_bowtie2_sort.log

      fi

      rm -rf ${s}_tmpdir

    # Calculate statistics from the alignments file

      echo $(date) 'Calculating statistics for '${s}

      java -Xmx3g -jar ${NGSEP} QualStats ${REF} ${s}_bowtie2_sorted.bam \
        >& ${s}_bowtie2_readpos.stats

      # java -Xmx3g -jar ${NGSEP} CoverageStats ${s}_bowtie2_sorted.bam \
      #   ${s}_bowtie2_coverage.stats >& ${s}_bowtie2_coverage.log

      echo $(date) ${s}' is DONE'

    done ) &

  done
  wait

  rm tmpList_*${pno}.tmp

  ################
  checkLastWord ${list[@]} _bowtie2_readpos.stats Bases
  ################

  echo -e '\n'$(date) 'Calculating Quality statistics for plate '${pno}

  cd ${WD}/mapping

  for col in {2..5}
  do

    touch tmpFile_${col}.tmp

    for sample in *_bowtie2_readpos.stats
    do

      head -n -3 ${sample} | cut -f${col} | paste tmpFile_${col}.tmp - > \
        all_pl${pno}_col${col}.txt
      cp all_pl${pno}_col${col}.txt tmpFile_${col}.tmp

    done

    sed -i "s/^[ \t]*//" all_pl${pno}_col${col}.txt
    rm tmpFile_${col}.tmp

    awk '{for(i=1;i<=NF;i++) t+=$i; print t; t=0}' all_pl${pno}_col${col}.txt > \
      sum_pl${pno}_col${col}.txt
    rm all_pl${pno}_col${col}.txt

    mv sum_pl${pno}_col${col}.txt ${WD}

  done

  cd ${WD}
  paste sum_pl${pno}_col2.txt sum_pl${pno}_col3.txt sum_pl${pno}_col4.txt \
    sum_pl${pno}_col5.txt > readPosStats_pl${pno}.tmp
  rm sum_pl*

  awk '{print $1/$3"\t"$2/$4}' readPosStats_pl${pno}.tmp > readPosStats_${pno}.txt

  echo -e '\nMapping on plate '${pno}' files seems to be completed\n'$(date)'\n'

  rm readPosStats_pl${pno}.tmp ${WD}/mapping/*bowtie2_readpos.stats
  printIt="stop"

else

  echo -e '\nMapping will not be executed this time\n'$(date)'\n'

fi


#### ------------------------------ ####
#### ----- Variant Discovery ------ ####


if [[ ${TASKS} == *V* ]]
then

  test ! -d ${WD}/mapping && \
    echo 'The mapping path '${WD}'/mapping does not exist !' && \
    exit 1

  cd ${WD}/mapping

  echo -e '\nStarting variant discovery on plate '${pno}' files\n'$(date)'\n'

  assignThreads ${list[@]} ${numThreads}

  for tmpFile in tmpList_*${pno}.tmp
  do

    # From the tmpList_XXX.tmp file (which is a chunk of the original list) 
    # put every line as an element of a bash array.
    myList=(`cat ${tmpFile} | tr '\n' ' '`)

    test ! ${printIt} && printThreads ${tmpFile}

  ( for s in ${myList[@]}
    do sleep 2

      echo $(date) 'Variant discovery in '${s}

      java -Xmx3g -jar ${NGSEP} FindVariants -h 0.0001 -maxBaseQS 30 \
        -minQuality 40 -maxAlnsPerStartPos 100 -ignore5 ${i5} -ignore3 ${i3} \
        -sampleId ${s} -knownSTRs ${STRs} ${REF} ${s}_bowtie2_sorted.bam \
        ${s}_bowtie2_NGSEP >& ${s}_bowtie2_NGSEP.log
      ${BGZIP} ${s}_bowtie2_NGSEP.vcf

    done ) &

  done
  wait

  rm tmpList_*${pno}.tmp

  ################
  checkLastWord ${list[@]} _bowtie2_NGSEP.log Completed
  ################

  echo -e '\nVariant-discovery on plate '${pno}' seems to be completed\n'$(date)'\n'

  printIt="stop"

else

  echo -e '\nThe FindVariants process will not be executed this time\n'$(date)'\n'

fi


#### ------------------------------ ####
#### ------- Quality-check -------- ####


if [[ ${TASKS} == *Q* ]]
then

  test ! -d ${WD}/mapping && \
    echo 'The mapping path '${WD}'/mapping does not exist !' && \
    exit 1

  cd ${WD}/mapping

  echo -e '\nStarting Quality-check on plate '${pno}'\n'$(date)'\n'

  # Avoid using too many threads here. This is a memory-demanding process
  assignThreads ${list[@]} 6

  for tmpFile in tmpList_*${pno}.tmp
  do

    # From the tmpList_XXX.tmp file (which is a chunk of the original list) 
    # put every line as an element of a bash array.
    myList=(`cat ${tmpFile} | tr '\n' ' '`)

    test ! ${printIt} && printThreads ${tmpFile}

  ( for s in ${myList[@]}
    do sleep 2

      echo $(date) 'Quality-check in '${s}

      java -Xmx20g -jar ${NGSEP} FindVariants -h 0.0001 -maxBaseQS 30 \
        -minQuality 60 -maxAlnsPerStartPos 100 -ignore5 ${i5} -ignore3 ${i3} \
        -sampleId ${s}-p${pno} -knownVariants ${myVariants} ${REF} \
        ${s}_bowtie2_sorted.bam ${s}-p${pno}_bowtie2_NGSEP_QC \
      >& ${s}-p${pno}_bowtie2_NGSEP_QC.log

      java -jar ${NGSEP} FilterVCF -q 60 -minI 1 ${s}-p${pno}_bowtie2_NGSEP_QC.vcf | \
      ${BGZIP} > ${s}-p${pno}_bowtie2_NGSEP_QC.vcf.gz

    done ) &

  done
  wait

  ################
  checkLastWord ${list[@]} -p${pno}_bowtie2_NGSEP_QC.log Completed
  ################

  echo -e '\nQuality-check on plate '${pno}' files seems to be completed\n'$(date)'\n'

  rm tmpList_*${pno}.tmp *bowtie2_NGSEP_QC.log *p${pno}_bowtie2_NGSEP_QC.vcf

else

  echo -e '\nThe Quality-check process will not be executed this time\n'$(date)'\n'

fi
