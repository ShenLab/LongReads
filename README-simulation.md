HapMut: Leveraging haplotype information in long reads for calling novel variants and mutations

The code base for HapMut has been cloned at the following location for convenience: /ifs/scratch/c2b2/ys_lab/yshen/HapMut/LongReads

This document attempts to explain how to use/modify the code base, locating datasets, running the tool etc.


There are two main parts of the code:

1) Simulation of datasets

2) Variant/mutation caller method

Simulation code is available under the directory: /ifs/scratch/c2b2/ys_lab/yshen/HapMut/LongReads while the methods
code base is under one of HMM, NB and Gibbs. Each directory contains code for the respective methods and apart from
that, the HMM directory also contains various scripts to plot/analyze the results.

The master pipeline is available at: /ifs/scratch/c2b2/ys_lab/yshen/HapMut/LongReads/pipeline.sh.
To run the pipeline for a single chromosome call:

  > nohup ./pipeline.sh 21 & (runs for chr21 only)
  
The main steps are:

  Simulate fastqs

    qsub -sync y -l mem=10G,time=4:: ./simulate_read.sh $snprate $errrate $coverage $readlen $region ${read_base} ${som_base}

  Map/Align reads
  
    qsub -sync y -l mem=8G,time=52:: ./long_read_map.sh ${read_base} $region
  
  Call preliminary variants
    
    qsub -sync y -l mem=1G,time=16:: ./snp_calling.sh ${read_base} ${vcf_base} $region
  
  Run respective method
    
    qsub -sync y -l mem=4G,time=32:: ./program.sh ${read_base} ${vcf_base} $region $iteration ${som_base} ${output_base} $pre $coverage $errrate
  
  Analyze results/Plots
    
    qsub -l mem=1G,time=2:: ./analysis.sh ${output_base}_${region} ${read_base}_${region} $coverage

Note that based on the variant caller method we adpot, the respective 'program' binary should be copied over to the
HMM directory before invoking the pipeline. The binary itself can be compiled in any of the directories by calling:
> make program

Also, in case the datasets are already simulated and variants called (mostly the case for testing), the first three
qsub steps can be simply commented out in the script.

Note the default parameters used. These can be appropriately modified.

snprate=0.02

errrate=0.02

coverage=10

readlen=2000

Analysis of results include plots of the mutation calling. When run as part of the pipeline, the plots are placed in
the HMM/output directory. Hence the plots for the above default parameters will be found at:

/ifs/scratch/c2b2/ys_lab/yshen/HapMut/LongReads/HMM/output/0.02_0.02_10_4000_21_som_ROC.png




Methods code base

Following are the main files in the methods codebase:

read.C: Contains data structure to store reads

snp.C:  Contains data structure to store snps

program.cpp: Contains code to read in the reads, snps, known snps etc. and call the respective variant caller method

hmm.C: Contains code to call haplotypes/genotypes/mutations and supporing routines

Any major changes to be made are expected to be made in hmm.C.
