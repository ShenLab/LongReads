#!/bin/sh
#$ -cwd
#$ -V
#$ -o recal.out
#$ -e recal.err
#$ -S /bin/bash
. /usr/prog/softenv/softenv.sh

JAVA="java -Xmx8000m -Djava.io.tmpdir=./tmp"
GATK="$JAVA -jar /usr/prog/onc/seqtools/GenomeAnalysisTK-2.3-9-ge5ebf34/GenomeAnalysisTK.jar"
GATK="$JAVA -jar /usr/prog/onc/seqtools/GenomeAnalysisTK-1.6-9-g47df7bb/GenomeAnalysisTK.jar"

PREFIX=$1
CHR="chr$2"

INP=${PREFIX}.realigned.bam
REF=reference/hg19.fasta
INDELVCF=reference/dbsnp137.vcf

$GATK \
  -T CountCovariates \
  -R $REF \
  -L $CHR \
  -knownSites $INDELVCF \
  -I $INP \
  -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate \
  -recalFile ${PREFIX}.recal_data.csv

$GATK \
  -T TableRecalibration \
  -L $CHR \
  -R $REF \
  -I $INP \
  -recalFile ${PREFIX}.recal_data.csv \
  -o ${PREFIX}.sorted.bam

mv ${PREFIX}.sorted.bai ${PREFIX}.sorted.bam.bai

##$GATK \
#    -T BaseRecalibrator \
#    -L $CHR \
 #   -I $INP \
#    -R $REF \
#    -knownSites $INDELVCF \
#    --maximum_cycle_value 5000 \
#    -o ${PREFIX}.recal_data.grp

##$GATK \
#    -T PrintReads \
#    -L $CHR \
#    -I $INP \
#    -R $REF \
#    -BQSR ${PREFIX}.recal_data.grp \
#    --disable_indel_quals \
#    -o ${PREFIX}.sorted.bam

