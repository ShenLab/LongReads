#!/bin/sh
#$ -cwd
#$ -V
#$ -o analysis.out
#$ -e analysis.err
#$ -S /bin/bash
. /usr/prog/softenv/softenv.sh

prefix=$1
inp=$2
coverage=$3

rm ${prefix}.a
qsub -sync y -l mem_token=1G -t 1-100 ./analyze.sh $prefix

#cp ../output/b.header ${prefix}.b
cp ../output/c.header ${prefix}.c

#sort -k4,4g ${prefix}.b[0-9]* >> ${prefix}.b
sort -k4,4g ${prefix}.c[0-9]* >> ${prefix}.c
#rm ${prefix}.b[0-9]*
rm ${prefix}.c[0-9]*

exit

hx=`grep HMM ${prefix}.a | grep pow | cut -d' ' -f 2`
hy=`grep HMM ${prefix}.a | grep pow | cut -d' ' -f 3`
sx=`grep Sam ${prefix}.a | grep pow | cut -d' ' -f 2`
sy=`grep Sam ${prefix}.a | grep pow | cut -d' ' -f 3`

# Code to obtain set of input read positions and their haplotype assignments
paste ${inp}.ibs ${prefix}.obs | awk 'BEGIN{stat=1;switch=0;known=0;ct=0;} {if($2!=$4) printf "reads %d and %d do not tally\n",$2,$4; else if($1!=$3) {printf "haplotypes %d and %d for reads %d and %d with %d known and %d snps do not tally\n", $1, $3, $2, $4, $6, $5; if(known>0) {if(stat!=2) switch++;ct++;} stat=2;} else {printf "haplotypes %d and %d for reads %d and %d with %d known and %d snps tally\n", $1, $3, $2, $4, $6, $5; if(known>0) {if(stat!=1) switch++;ct++;} stat=1; sum++;}known=$6; } END{printf "Total %d tally from %d: %f with %d switch errors of %d and switch error rate %f\n", sum, NR, (sum/NR)*100, switch, ct, (switch/ct)*100;}' >> ${prefix}.a

cp base.gnu ${prefix}.gnu
echo "set xlabel '(1 - Specificity)'" >> ${prefix}.gnu
echo "set ylabel 'Sensitivity'" >> ${prefix}.gnu
echo "set output '${prefix}_som_ROC.png'" >> ${prefix}.gnu
echo "plot '${prefix}.c' using 2:3 w l ls 1, \\" >> ${prefix}.gnu
echo "	'${prefix}.c' using 6:7 w l ls 2" >> ${prefix}.gnu
~/bin/gnuplot ${prefix}.gnu

cp base.gnu ${prefix}.gnu
echo "set xlabel '(1 - Accuracy)'" >> ${prefix}.gnu
echo "set ylabel 'Power'" >> ${prefix}.gnu
echo "set output '${prefix}_som_Power.png'" >> ${prefix}.gnu
echo "plot '${prefix}.c' using 1:3 w l ls 1, \\" >> ${prefix}.gnu
echo "	'${prefix}.c' using 5:7 w l ls 2" >> ${prefix}.gnu
~/bin/gnuplot ${prefix}.gnu

cp base.gnu ${prefix}.gnu
echo "set xlabel '(1 - Full Accuracy)'" >> ${prefix}.gnu
echo "set ylabel 'Power'" >> ${prefix}.gnu
echo "set output '${prefix}_som_FDR.png'" >> ${prefix}.gnu
echo "plot '${prefix}.c' using 8:3 w l ls 1, \\" >> ${prefix}.gnu
echo "	'${prefix}.c' using 9:7 w l ls 2" >> ${prefix}.gnu
~/bin/gnuplot ${prefix}.gnu

cp base.gnu ${prefix}.gnu
echo "set xlabel '(1 - Normalized Specificity)'" >> ${prefix}.gnu
echo "set ylabel 'Sensitivity'" >> ${prefix}.gnu
echo "set output '${prefix}_som_Norm_ROC.png'" >> ${prefix}.gnu
echo "plot '${prefix}.c' using 10:3 w l ls 1, \\" >> ${prefix}.gnu
echo "	'${prefix}.c' using 11:7 w l ls 2" >> ${prefix}.gnu
~/bin/gnuplot ${prefix}.gnu

# Code to find point of maximum distance on the plot from a (random) diagonal
## echo "Farthest point from diagonal" >> ${prefix}.a
#awk 'function ddist(x,y) {min=100;x=x*100;y=y*100; for(i=1;i<=100;i++) {if(sqrt((x-i)*(x-i)+(y-i)*(y-i))<min) min=sqrt((x-i)*(x-i)+(y-i)*(y-i));} return min;} BEGIN {hpd=0;hpx=0;hpy=0;hpi=0;hrd=0;hrx=0;hry=0;hri=0;spd=0;spx=0;spy=0;spi=0;srd=0;srx=0;sry=0;sri=0;} {if(ddist($1,$3)>hpd) {hpd=ddist($1,$3);hpx=$1;hpy=$3;hpi=$4;} if(ddist($2,$3)>hrd) {hrd=ddist($2,$3);hrx=$2;hry=$3;hri=$4;} if(ddist($5,$7)>spd) {spd=ddist($5,$7);spx=$5;spy=$7;spi=$4;} if(ddist($6,$7)>srd) {srd=ddist($6,$7);srx=$6;srx=$7;sri=$4;}} END{print "HMM(pow,acc)",hpx,hpy,hpd,hpi;print "Sam(pow,acc)",spx,spy,spd,spi;print "HMM(spc,sen)",hrx,hry,hrd,hri;print "Sam(spc,sen)",srx,sry,srd,sri;}' ${prefix}.b >> ${prefix}.a
## echo >> ${prefix}.a

# Maybe I'll combine these in the super master pipeline script to plot for all simulation parameters in single plot

#cp base.gnu ${prefix}.gnu
#echo "set xlabel '(1 - Specificity)'" >> ${prefix}.gnu
#echo "set ylabel 'Sensitivity'" >> ${prefix}.gnu
#echo "set output '${prefix}_ROC.png'" >> ${prefix}.gnu
#echo "plot '${prefix}.b' using 2:3 w l ls 1, \\" >> ${prefix}.gnu
#echo "	'${prefix}.b' using 6:7 w l ls 2" >> ${prefix}.gnu
# echo "	'${prefix}.c' using 2:3 w l ls 3" >> ${prefix}.gnu
#gnuplot ${prefix}.gnu

#cp base.gnu ${prefix}.gnu
#echo "set xlabel '(1 - Accuracy)'" >> ${prefix}.gnu
#echo "set ylabel 'Power'" >> ${prefix}.gnu
#echo "set output '${prefix}_Power.png'" >> ${prefix}.gnu
#echo "plot '${prefix}.b' using 1:3 w l ls 1, \\" >> ${prefix}.gnu
#echo "	'${prefix}.b' using 5:7 w l ls 2, \\" >> ${prefix}.gnu
#echo "	\"<echo '$hx $hy'\" w p, \\" >> ${prefix}.gnu
#echo "	\"<echo '$sx $sy'\" w p" >> ${prefix}.gnu
#gnuplot ${prefix}.gnu

# cd output
# sr=`basename $prefix | cut -d'_' -f1`
# ir=`basename $prefix | cut -d'_' -f2`
# cov=`basename $prefix | cut -d'_' -f3`
# rl=`basename $prefix | cut -d'_' -f4`

# prev=`echo "scale=2; (100-2*$coverage)/100" | bc -l`
# Rscript plot_novel.R $sr $ir $cov $rl
# Rscript plot_som.R $sr $ir $cov $rl $prev
#  ./som.sh $coverage
# cd ..

# Code to find sensitivity to accuracy ratio based on maximum value of genotype - single point in plot
# Maybe split for known and novel hets later?
# echo "Sensitivity and accuracy (1-specificity) based on highest genotype value" > ${prefix}.a
#awk '{if($2==1||$2==2) {tot++;if($3==1) sum++;} else {pot++; if($3!=1) pum++;}} END{sens = sum/tot; spc = pum/pot; print 1-sens, "\t", spc;}' ${prefix}.sta >> ${prefix}.a
#echo >> ${prefix}.a

# Code to find area under curve of sensitivity to accuracy plot
# calculated using G=2AUC-1 where, G = 1 - sum(X(k)-X(k-1)*Y(k)+Y(k-1)), k= {2,..,n}
#echo "Area under curve" >> ${prefix}.a
#awk '{a[NR]=$1;b[NR]=$2;} END{for(i=NR-1;i>=0;i--) {sum += ((a[i]-a[i+1])*(b[i]+b[i+1]))} print 1-sum/2;}' ${prefix}.b >> ${prefix}.a
#echo >> ${prefix}.a

