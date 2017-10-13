#!/bin/bash

TASSEL=/home/dariza/software/TASSEL5/run_pipeline.pl
HMP=/bioinfo1/projects/bean/VEF/GWAS/noMeso_mlm-PCA/VEF_noMeso_annotated_repMasked_q40_s_fi_maf05_oh06_i210_imputed_hmp.txt
PHEN=$1
OUT=$2
LOG=$3

# ${TASSEL} -h ${HMP} -sortPositions -KinshipPlugin -method Centered_IBS -endPlugin -export VEF_annotated_repMasked_q40_s_fi_maf05_oh06_i290_kinship.txt -exportType SqrMatrix

Kinship=/bioinfo1/projects/bean/VEF/GWAS/noMeso_mlm-PCA/VEF_noMeso_annotated_repMasked_q40_s_fi_maf05_oh06_i210_imputed_kinship.txt

# ${TASSEL} -fork1 -h ${HMP} -sortPositions -PrincipalComponentsPlugin -covariance true -ncomponents 10 -reportEigenvalues true -reportEigenvectors true -endPlugin -export VEF_annotated_repMasked_q40_s_fi_maf05_oh06_i290_PCA.txt -runfork1

PCA=/bioinfo1/projects/bean/VEF/GWAS/noMeso_mlm-PCA/VEF_noMeso_annotated_repMasked_q40_s_fi_maf05_oh06_i210_imputed_PCA1.txt

${TASSEL} -Xms512m -Xmx10g -log ${LOG} -fork1 -h ${HMP} -sortPositions -fork2 -r ${PHEN} -fork3 -q ${PCA} -fork4 -k ${Kinship} -combine5 -input1 -input2 -input3 -intersect -combine6 -input5 -input4 -mlm -mlmVarCompEst P3D -mlmCompressionLevel None -export ${OUT} -runfork1 -runfork2 -runfork3 -runfork4
