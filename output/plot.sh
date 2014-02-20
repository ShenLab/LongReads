#!/bin/sh
#$ -cwd

for prev in 10.0 25.0
do
  for rlen in 500 1000 2000 4000
  do
    for cov in 10 20 30 50
    do
      echo $prev $rlen $cov
      Rscript plot_ROC.R $prev $rlen $cov
    done
  done
done

for prev in 10.0 25.0
do
  for errate in 1.000 2.00 4.00
  do
    grep "$prev%" table.tsv | grep "$errate%" | grep 500 | sort -k3,3g > _500_.txt
    grep "$prev%" table.tsv | grep "$errate%" | grep 1000 | sort -k3,3g > _1000_.txt
    grep "$prev%" table.tsv | grep "$errate%" | grep 2000 | sort -k3,3g > _2000_.txt
    grep "$prev%" table.tsv | grep "$errate%" | grep 4000 | sort -k3,3g > _4000_.txt
    echo "jpeg(filename='TP_${prev}_${errate}.jpeg')" > plot_TP.R
    echo "jpeg(filename='FP_${prev}_${errate}.jpeg')" > plot_FP.R
    echo "ttl = \"err=${errate}%, prev=${prev}%"\" >> plot_TP.R
    echo "ttl = \"err=${errate}%, prev=${prev}%"\" >> plot_FP.R
    cat base_TP.R >> plot_TP.R
    cat base_FP.R >> plot_FP.R
    Rscript plot_TP.R
    Rscript plot_FP.R
  done
done

