#!/bin/bash

list=$1

WD=/bioinfo2/projects/bean_tmp/organizing_data_tmp/ADP
plate=ADP

cd ${WD}/mapping

for col in {2..5}
do

  touch tmpFile_${col}.tmp

  for sample in *_bowtie2_readpos.stats
  do

    head -n -3 ${sample} | cut -f${col} | paste tmpFile_${col}.tmp - > all_pl${plate}_col${col}.txt
    cp all_pl${plate}_col${col}.txt tmpFile_${col}.tmp

  done

  sed -i "s/^[ \t]*//" all_pl${plate}_col${col}.txt
  rm tmpFile_${col}.tmp

  awk '{for(i=1;i<=NF;i++) t+=$i; print t; t=0}' all_pl${plate}_col${col}.txt > sum_pl${plate}_col${col}.txt
  rm all_pl${plate}_col${col}.txt

  mv sum_pl${plate}_col${col}.txt ${WD}


done

cd ${WD}
paste sum_pl${plate}_col2.txt sum_pl${plate}_col3.txt sum_pl${plate}_col4.txt sum_pl${plate}_col5.txt > readPosStats_pl${plate}.tmp
rm sum_pl*

awk '{print $1/$3"\t"$2/$4}' readPosStats_pl${plate}.tmp > readPosStats_pl${plate}.txt
rm readPosStats_pl${plate}.tmp


