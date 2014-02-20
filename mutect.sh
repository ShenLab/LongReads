#!/bin/sh
#$ -cwd
#$ -V
#$ -o file.out
#$ -e file.err
#$ -S /bin/bash
. /usr/prog/softenv/softenv.sh

JAVA="java -Xmx8000m -Djava.io.tmpdir=./tmp"
GATK="$JAVA -jar /usr/prog/onc/seqtools/GenomeAnalysisTK-2.3-9-ge5ebf34/GenomeAnalysisTK.jar"
GATK="$JAVA -jar /usr/prog/onc/seqtools/GenomeAnalysisTK-1.6-9-g47df7bb/GenomeAnalysisTK.jar"
MUTECT="$JAVA -jar mutect/muTect-1.1.4.jar"

PREFIX=reads/$1_$1_$2_$3_$5_$4
OUTPRE=output/$1_$1_$2_$3_$5_$4
CHR=chr$4

INP=${PREFIX}.sorted.bam
REF=reference/hg19.fasta
INDELVCF=reference/dbsnp137.vcf

$MUTECT \
 --analysis_type MuTect \
 --reference_sequence $REF \
 --dbsnp $INDELVCF \
 --cosmic mutect/hg19_cosmic_v54_120711.vcf \
 --intervals $CHR \
 --input_file:tumor $INP \
 --out ${OUTPRE}.call_stats.txt \
 --coverage_file ${OUTPRE}.coverage.wig.txt

som=`cat snps/som_$1_$1_$2_$3_$5_$4.list`
germ=`cut -f 2 snps/snp_$4.list`
mutkeep=`grep KEEP ${OUTPRE}.call_stats.txt | cut -f 2`
mutall=`cut -f 2 ${OUTPRE}.call_stats.txt`
somall=`echo $mutall $som | tr ' ' '\n' | sort -g | uniq -c | awk '$1==2' | wc -l`

mutkeepgerm=`echo $germ $mutkeep | tr ' ' '\n' | sort -g | uniq -c | awk '$1==2' | wc -l`
mutallgerm=`echo $germ $mutall | tr ' ' '\n' | sort -g | uniq -c | awk '$1==2' | wc -l`

mutkeepsom=`echo $som $mutkeep | tr ' ' '\n' | sort -g | uniq -c | awk '$1==2' | wc -l`
mutallsom=`echo $som $mutall | tr ' ' '\n' | sort -g | uniq -c | awk '$1==2' | wc -l`

ktp=`echo $mutkeep $som | tr ' ' '\n' | sort -g | uniq -c | awk '$1==2' | wc -l`
ktot=`echo $mutkeep | tr ' ' '\n' | wc -l`
let "kfp=${ktot}-$ktp-$mutkeepgerm"

atp=`echo $mutall $som | tr ' ' '\n' | sort -g | uniq -c | awk '$1==2' | wc -l`
atot=`echo $mutall | tr ' ' '\n' | wc -l`
let "afp=${atot}-$atp-$mutallgerm"

echo "" > ${OUTPRE}.m
echo "Keep $ktp $kfp" >> ${OUTPRE}.m
echo "All $atp $afp" >> ${OUTPRE}.m
somall=160
echo "Num Sensitivity FalsePositives TruePositives" | tr ' ' '\t' > ${OUTPRE}.n
for i in `seq 1 100`
do
	mutrockeep=`grep KEEP ${OUTPRE}.call_stats.txt | cut -f 2,17 | awk -v thresh=$i '$2>thresh' | cut -f 1`
	mutrockeepgerm=`echo $germ $mutrockeep | tr ' ' '\n' | sort -g | uniq -c | awk '$1==2' | wc -l`
	rktp=`echo $mutrockeep $som | tr ' ' '\n' | sort -g | uniq -c | awk '$1==2' | wc -l`
	sktp=`echo "scale=2; $rktp/$somall" | bc`
	rktot=`echo $mutrockeep | tr ' ' '\n' | wc -l`
	let "rkfp=${rktot}-$rktp-$mutrockeepgerm"
	num=`echo "scale=2; $i/100" | bc`
	echo "$num $sktp $rkfp $rktp" | tr ' ' '\t' >> ${OUTPRE}.n
done

