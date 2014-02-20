#!/bin/sh
#$ -cwd
#$ -V
#$ -o file.out
#$ -e file.err
#$ -S /bin/bash
. /usr/prog/softenv/softenv.sh

snprate=0.02
errrate=0.02
coverage=50
readlen=4000
region=21
freq=15
iteration=1

snprate=$1
errrate=$2
coverage=$3
readlen=$4
region=$5
freq=$6

base=$PWD
pre=${snprate}_${errrate}_${coverage}_${readlen}_$freq
read_base=${base}/reads/${pre}
vcf_base=${base}/vcfs/${pre}
som_base=${base}/snps/som_${pre}

#gcc -g simulate.c -o simulate 

if [[ $region -ne 0 ]]
then
	qsub -sync y -l mem_token=10G ./simulate_read.sh $snprate $errrate $coverage $readlen $region ${read_base} ${som_base} $freq
	echo simulation complete

	qsub -sync y -l mem_token=20G ./long_read_map.sh ${read_base} $region
	echo alignment complete

	qsub -sync y -l mem_token=1G ./snp_calling.sh ${read_base} ${vcf_base} $region $snprate $freq
	echo snp calling complete

	qsub -l mem_token=9G ./mutect.sh $snprate $coverage $readlen $region $freq
	echo mutect kicked in

	cd Gibbs 
	output_base=${base}/output/${pre}
	qsub -sync y -l mem_token=4G ./program.sh ${read_base} ${vcf_base} $region $iteration ${som_base} ${output_base} $pre $coverage $errrate $freq
	echo HMM complete
	qsub -l mem_token=1G ./analysis.sh ${output_base}_${region} ${read_base}_${region} $coverage $freq
else
	qsub -sync y -t 1-22 -l mem_token=10G ./simulate_read.sh $snprate $errrate $coverage $readlen $region ${read_base} ${som_base}
	echo simulation complete

	qsub -sync y -t 1-22 -l mem_token=8G ./long_read_map.sh ${read_base} $region

	rm ${read_base}.ibs
	for i in `seq 1 22`
	do
		echo ${read_base}_${i}.ibs
		cat ${read_base}_${i}.ibs >> ${read_base}.ibs
		rm ${read_base}_${i}.ibs
	done
	echo alignment complete

	qsub -sync y -t 1-22 -l mem_token=1G ./snp_calling.sh ${read_base} ${vcf_base} $region
	echo snp calling complete

	cd Gibbs
	output_base=${base}/output/${pre}_${iteration}
	qsub -sync y -t 1-22 -l mem_token=4G ./program.sh ${read_base} ${vcf_base} $region $iteration ${som_base} ${output_base} $pre $coverage $errrate

	rm ${output_base}.obs
	rm ${output_base}.sta
	for i in `seq 1 22`
	do
		echo $i
		cat ${output_base}_${i}.obs >> ${output_base}.obs
		cat ${output_base}_${i}.sta >> ${output_base}.sta
		rm ${output_base}_${i}.obs
		rm ${output_base}_${i}.sta
	done
	echo HMM complete

	qsub -l mem_token=1G ./analysis.sh ${output_base} ${read_base} $coverage
fi

cd ..

echo Analysis complete

