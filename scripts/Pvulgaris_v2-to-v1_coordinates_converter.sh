#!/bin/bash

chromosome=$1; # Specify the chromosome ID. Use the same ID as it appears on the fasta reference file.
start=$2; # Start position to be converted
end=$3; # End position of the range to be converted. If you only want to convert a single position (instead of a range), leave the 3rd argument empty or use the same value as in 'start'
# fastafile=$4; 
numThreads=10

blast=/data/software/ncbi-blast-2.2.27+/bin
samtools=/data/software/samtools/samtools-1.2/samtools


REFv1=/data/references/bean/v1.0/Pvulgaris218_phytozome.fa; # This file must have a '.fai' fasta index in the same path
REFv2=/data/references/bean/v2.1/blast/Pvulgaris_442_v2.0.fa; # This file must be a BLAST indexed database. Check BLAST's makeblastdb

if [[ ! -z "$2" ]]
then 

  if [[ ${start} == ${end} ]] || [[ -z "$3" ]]; then end=`expr ${start} + 500`; fi

  ${samtools} faidx ${REFv1} ${chromosome}:${start}-${end} > ${chromosome}-${start}-${end}.fasta
  fastafile=${chromosome}-${start}-${end}.fasta

else

  fastafile=$1

fi

${blast}/blastn -query ${fastafile} -db ${REFv2} -out ${chromosome}-${start}-${end}.outfmt6 -evalue 1e-10 -num_threads ${numThreads} -outfmt 0 -task megablast -best_hit_score_edge 0.05 -best_hit_overhang 0.25 -max_target_seqs 5

cat ${chromosome}-${start}-${end}.outfmt6

rm ${chromosome}-${start}-${end}.fasta ${chromosome}-${start}-${end}.outfmt6
