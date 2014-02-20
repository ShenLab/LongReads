#!/bin/sh
#$ -cwd

suf1=$1
suf2=$2

for pre in `cat param.list`
do
	pre="bkup/$pre"
	i=`basename $pre`
	echo "${i}	True Positives	False positives"
	tp=`awk '$13==1&&$3==0' ${pre}.sta.$suf1 | wc -l`
	fp=`awk '$13==0&&$3==0' ${pre}.sta.$suf1 | wc -l`
	echo "Total	$tp	$fp"
	mktp=`grep Keep ${pre}.m.$suf2 | cut -d' ' -f 2`
	mkfp=`grep Keep ${pre}.m.$suf2 | cut -d' ' -f 3`
	echo "Mutect-Keep	$mktp	$mkfp"
	matp=`grep All ${pre}.m.$suf2 | cut -d' ' -f 2`
	mafp=`grep All ${pre}.m.$suf2 | cut -d' ' -f 3`
	echo "Mutect-All	$matp	$mafp"
	tpg=`awk '$13==1&&$3==0&&$14==1' ${pre}.sta.$suf1 | wc -l`
	fpg=`awk '$13==0&&$3==0&&$14==1' ${pre}.sta.$suf1 | wc -l`
	echo "==1	$tpg	$fpg"
	tpg=`awk '$13==1&&$3==0&&$14>.99' ${pre}.sta.$suf1 | wc -l`
	fpg=`awk '$13==0&&$3==0&&$14>.99' ${pre}.sta.$suf1 | wc -l`
	echo ">.99	$tpg	$fpg"
	tpg=`awk '$13==1&&$3==0&&$14>.75' ${pre}.sta.$suf1 | wc -l`
	fpg=`awk '$13==0&&$3==0&&$14>.75' ${pre}.sta.$suf1 | wc -l`
	echo ">.75	$tpg	$fpg"
	tpg=`awk '$13==1&&$3==0&&$14>.5' ${pre}.sta.$suf1 | wc -l`
	fpg=`awk '$13==0&&$3==0&&$14>.5' ${pre}.sta.$suf1 | wc -l`
	echo ">.5	$tpg	$fpg"
	tpg=`awk '$13==1&&$3==0&&$14>.25' ${pre}.sta.$suf1 | wc -l`
	fpg=`awk '$13==0&&$3==0&&$14>.25' ${pre}.sta.$suf1 | wc -l`
	echo ">.25	$tpg	$fpg"
	tpg=`awk '$13==1&&$3==0&&$14>.1' ${pre}.sta.$suf1 | wc -l`
	fpg=`awk '$13==0&&$3==0&&$14>.1' ${pre}.sta.$suf1 | wc -l`
	echo ">.1	$tpg	$fpg"
	tpg=`awk '$13==1&&$3==0&&$14>.01' ${pre}.sta.$suf1 | wc -l`
	fpg=`awk '$13==0&&$3==0&&$14>.01' ${pre}.sta.$suf1 | wc -l`
	echo ">.01	$tpg	$fpg"
	tpg=`awk '$13==1&&$3==0&&$14>0' ${pre}.sta.$suf1 | wc -l`
	fpg=`awk '$13==0&&$3==0&&$14>0' ${pre}.sta.$suf1 | wc -l`
	echo ">0	$tpg	$fpg"
done

