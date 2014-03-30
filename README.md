HapMut: Leveraging haplotype information in long reads for calling somatic mutations
====================================================================================


Overview: 
---------

HapMut is a tool used to identify somatic mutations from tumor samples, without matched normals. HapMut works best on bam files with long reads (500bp+) generated from a second-generation (e.g., Illumina HMiSeqTM ) or third-generation (e.g., nanopore technologies) sequencing technologies. Long reads of length 500bp and more harbor multiple polymorphic sites with reasonable frequency, and contain haplotype information to infer a true mutation from an error. HapMut uses a Naive Bayes Classifier to leverage such haplotype information and call the haplotypes in a sample. Here we assume a cancer sequencing study is designed as following: individual cancer cells are separated from a single tissue or multiple tissues; each cell is sequenced using a third-generation technology without pre-amplification or a second-generation technology with pre-amplification.

HapMut takes as input a VCF file containing a preliminary set of variant calls, and refines the genotypes of the variants to infer true somatic mutations. The preliminary set of variants maybe generated using tools such as samtools or GATK. HapMut also requires a [dbSNP file](ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/) as input. The dbSNP calls serve two purposes: 1) It provides a set of "known" variants or polymorphic sites with high confidence for haplotype calling. 2) Since the dbSNP calls are more likely to be germline than somatic, these calls are filtered out from the preliminary set in the context of unmatched normals. Filtering out such calls is very useful in the absence of a matched normal sample since otherwise all germline calls would also be reported as somatic. HapMut also employs other filters that distinguish sequencing or alignment artifacts from true mutations.


Source Package:
---------------

The source code for HapMut is freely available for download at https://github.com/joshuashen/HapMut. The code can be downloaded using the following git command:

git clone https://github.com/joshuashen/HapMut

The code is written in C++ and tested to work on a Unix platform.


Compilation:
------------

HapMut requires g++, the standard C/C++ libraries and also the [Bamtools API library](https://github.com/pezmaster31/bamtools) for compilation. A makefile is provided to compile HapMut using the necessary libraries. Before proceeding with compilation, the path to the header files and the shared libraries of Bamtools API has to be specified in the makefile. HapMut can then be compiled using the 'make hapmut' command in the same directory as the source files. The generated binary can then be placed in your exported PATH location or any other local folder.


Usage:
------

HapMut requires four minimum parameters to run: 1) bam file from a single sample. 2) VCF file containing preliminary variant calls. 3) dbsNP VCF file. 4) Output base path. Other parameters can be provided such as the chromosome or region of interest, maximum read length and expected error rate of reads. For all hapmut parameters and their description, simply run ./hapmut after compilation.


Mutation calling using HapMut:
------------------------------

Given a fastq file to begin with, somatic mutations can be called by HapMut using the following steps:

 1. **Read alignment** - First align the reads to the reference genome using a suitable aligner such as [bwa](http://bio-bwa.sourceforge.net/bwa.shtml) or [Stampy](http://www.well.ox.ac.uk/project-stampy). If using bwa to align the reads, it is recommended to perform indel realignment and base quality recalibration using [GATK](www.broadinstitute.org/gatk/) as mentioned in the [GATK best practices document](http://www.broadinstitute.org/gatk/guide/best-practices). This will give you a bam file containing the aligned reads.
 2. **Sorting and duplicates** - Sort the bam file using [samtools](http://samtools.sourceforge.net/) or another tool and index the sorted bam file. If you used PCR to prepare your sequencing library, mark the duplicates in the bam file using a tools such as [Picard](http://picard.sourceforge.net/command-line-overview.shtml#MarkDuplicates) and index the new bam file as well.
 3. **Preliminary variant calling** - Before using HapMut, obtain a preliminary set of variant calls using a variant calling tool such as samtools or GATK. It is recommended NOT to keep any stringent parameters when calling variants, or perform any variant filtration post the variant calling. If using bcftools view after samtools mpileup to prepare the input VCF, the default p-value of 0.5 can be increased to 1.0 using the '-p' switch to consider an exhaustive set of potential mutations.
 4. **HapMut mutation calling** - We are now ready to call HapMut to refine the genotypes and report somatic mutation calls. HapMut can be run as described in the previous step.


Output mutations:
-----------------

HapMut outputs the list of input variant calls with the probability of each variant being somatic. Based on a certain threshold, HapMut determines whether a variant call is somatic or not. Candidate sites also belonging to dbSNP are not output as they are considered germline.


Contact:
--------

For any questions or issues with HapMut, please email to ys2411@columbia.edu.

