#!/bin/sh
#$ -cwd

pre=$1

# awk '$13==1&&$14>$9' $file | wc -l
# awk '$13==1&&$14==$9' $file | wc -l
# awk '$13==1&&$14<$9' $file | wc -l
# awk '$13==0&&$14<$9' $file | wc -l
# awk '$13==0&&$14==$9' $file | wc -l
# awk '$13==0&&$14>$9' $file | wc -l

#awk '{if($13==1) {if($14>.99) f99++; if($14>.75) f75++; if($14>.5) f50++; if($14>.25) f25++; if($14>.1) f10++; if($14>.01) f01++; f++;} else if($13==0) {if($14>.99) s99++; if($14>.75) s75++; if($14>.5) s50++; if($14>.25) s25++; if($14>.1) s10++; if($14>.01) s01++; s++;} } END{ print "Total 99 75 50 25 10 1";print f,f99,f75,f50,f25,f10,f01; print s,s99,s75,s50,s25,s10,s01;}' ${pre}.sta

echo "Read_length Error_rate Coverage Allelic_Freq HapMut_Sensitivity Hapmut_FP Mutect_Sensitivity Mutect_FP" | tr ' ' '\t'

for pre in `ls 17/0*15_21.sta | cut -d'.' -f 1-3`
do
hm=`awk '{if($13==1&&$3==0) {if($14>.5) f50++; f++;} else if($13==0&&$3==0) {if($14>.5) s50++;s++;} } END{ printf "%d %.1f %d\n", f,f50*100/160,s50;}' ${pre}.sta`
hs=`echo $hm | cut -d' ' -f 2`
hf=`echo $hm | cut -d' ' -f 3`

#total=`echo $hm | cut -d' ' -f 1`
total=160

mt=`grep Keep ${pre}.m`
mp=`echo $mt | cut -d' ' -f 2`
mf=`echo $mt | cut -d' ' -f 3`
ms=`echo "scale=1; ${mp}*100/$total" | bc`

er=`echo $pre | cut -d'_' -f 1`
err=`echo "scale=1; $er*200" | bc`
len=`echo $pre | cut -d'_' -f 4`
cov=`echo $pre | cut -d'_' -f 3`

echo $len $err% $cov 10.0% $hs $hf $ms $mf | tr ' ' '\t'
done

for pre in `ls 17/0*35_21.sta | cut -d'.' -f 1-3`
do
hm=`awk '{if($13==1&&$3==0) {if($14>.5) f50++; f++;} else if($13==0&&$3==0) {if($14>.5) s50++;s++;} } END{ printf "%d %.1f %d\n", f,f50*100/160,s50;}' ${pre}.sta`
hs=`echo $hm | cut -d' ' -f 2`
hf=`echo $hm | cut -d' ' -f 3`

#total=`echo $hm | cut -d' ' -f 1`
total=160

mt=`grep Keep ${pre}.m`
mp=`echo $mt | cut -d' ' -f 2`
mf=`echo $mt | cut -d' ' -f 3`
ms=`echo "scale=1; ${mp}*100/$total" | bc`

er=`echo $pre | cut -d'_' -f 1`
err=`echo "scale=1; $er*200" | bc`
len=`echo $pre | cut -d'_' -f 4`
cov=`echo $pre | cut -d'_' -f 3`

echo "$len $err% $cov 25.0% $hs $hf $ms $mf" | tr ' ' '\t'
done

