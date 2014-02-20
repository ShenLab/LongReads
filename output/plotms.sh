#!/bin/sh
#$ -cwd

for prev in 10.0 25.0
do
  for cov in 10 20 30 50
  do
echo $prev $rlen $cov
    Rscript plotms_ROC.R $prev $cov
  done
done

Rscript plotms_TP.R
Rscript plotms_FP.R

