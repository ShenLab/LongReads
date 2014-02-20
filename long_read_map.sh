#!/bin/sh
#$ -cwd
#$ -V
#$ -o long.out
#$ -e long.err
#$ -S /bin/bash
. /usr/prog/softenv/softenv.sh

base=$1
region=$2

if [[ $region -eq 0 ]]
then
	region=$SGE_TASK_ID
fi

prefix=${base}_${region}

PICARD=/usr/prog/picard-tools/1.52/
JAVA="java -Xmx19000m -Djava.io.tmpdir=./tmp"

# One-time calling
# bwa index -a bwtsw bcm_hg18.fasta

# Defaults are: -b 3 -q 5 -r 2
# Expt range: -b [4-6] -q [2-4] -r [1-2]

#bwa bwasw -b 4 -q 4 -r 4 -w 500 reference/hg19.fasta ${prefix}.fq  > ${prefix}_temp.samd
#bwa mem -H -t 4 -w 50 -O 4 -E 6 -R "@RG\tID:hmID\tSM:hmSM\tLB=hmLB\tPL=illumina\tPU=hmPU" reference/hg19.fasta ${prefix}.fq  > ${prefix}_temp.samd

bwa mem -H -w 50 -v 0 reference/hg19.fasta ${prefix}.fq  > ${prefix}_temp.samd

#awk -v reg="chr${region}" 'function rem(cigar) {num=0;for(it=1;it<=length(cigar);it++) {sb=substr(cigar,it,1); if(sb ~ /[0-9]/) {num=10*num+sb;} else { if(sb!="H") return 0; else return num;}}} {if((NF==3)||!match($10,"N")&&$4>0&&$5>=10&&$3==reg&&$6!="*"&&!match($6,"H")&&!match($6,"S")) print}' ${prefix}_temp.samd > ${prefix}_temp.sam
awk -v reg="chr${region}" 'function rem(cigar) {num=0;for(it=1;it<=length(cigar);it++) {sb=substr(cigar,it,1); if(sb ~ /[0-9]/) {num=10*num+sb;} else { if(sb!="H") return 0; else return num;}}} {if((NF==3)||!match($10,"N")&&$4>0&&$3==reg&&$6!="*") print}' ${prefix}_temp.samd > ${prefix}_temp.sam
samtools calmd -S ${prefix}_temp.sam reference/hg19.fasta > ${prefix}.sam

rm ${prefix}_temp.samd ${prefix}_temp.sam
samtools view -bS ${prefix}.sam | samtools sort -m 2000000000 - ${prefix}.reheader
samtools reheader reads/header.sam ${prefix}.reheader.bam > ${prefix}.sort.bam
samtools index ${prefix}.sort.bam

$JAVA -jar ${PICARD}/AddOrReplaceReadGroups.jar I= ${prefix}.sort.bam O= ${prefix}.rg.bam SORT_ORDER=coordinate RGID=HapMut1 RGLB=human RGPL=Illumina RGSM=HMSM1 RGPU=HMPU1
samtools index ${prefix}.rg.bam

qsub -sync y gatk_realign.sh ${prefix} $region
qsub -sync y gatk_recalibrate.sh ${prefix} $region

#cp ${prefix}.realigned.bam ${prefix}.sorted.bam
#samtools index ${prefix}.sorted.bam

samtools view -h ${prefix}.sorted.bam > ${prefix}.sorted.sam
grep -v ^@ ${prefix}.sorted.sam | cut -f 1,4 | cut -d ':' -f 3 > ${prefix}.ibs

#rm ${prefix}.fq
rm ${prefix}.reheader.bam ${prefix}.forRealigner.intervals ${prefix}.recal_data.csv ${prefix}.realigned.bai ${prefix}.realigned.bam ${prefix}.rg.bam ${prefix}.rg.bam.bai ${prefix}.sort.bam ${prefix}.sort.bam.bai ${prefix}.sam

