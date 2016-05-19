# *STR-RNA-MLE*

A program written in R that can compute the maximum likelihood estimate of the rates of RNA-DNA differences and reverse transcription errors (and associated expansion or contraction probabilities) from short tandem repeat profiles generated with RNA-seq.


## Overview
The scripts to perform full MLE and lumping MLE are designed to estimate RNA-DNA difference (RDD) and reverse transcriptase error rates/patterns of short tandem repeats (STRs) simultaneously using genome and transcriptome data. These algorithms estimate the following parameters:
- RDD rate = The proportion of the STR repeat number in RNA that are deviated from the STR repeat number of template DNA due to either RNA editing or transcription slippage.
- RDD expansion probability = The probability that the RDD are expansion.
- RT error rate = The proportion of the STR repeat number in cDNA that are deviated from the STR repeat number of RNA due to reverse transcription process.
- RT expansion probability = The probability that the RT errors are expansion.
 
These algorithms will calculate the probability that the STR length profiles of all loci that share the same genotype were generated from specific set of RDD rate, RDD expansion probability, RT error rate, RT error expansion probability. Each set of these four parameters is chosen by L-BFGS-B algorithm from optim function in R. At the end of the process, the algorithm will report the optimal set of parameters that maximize the likelihood that the data were generated. The sequencing error rates from Fungtammasan et al. 2015 were used as constant rates of transition from cDNA to sequencing read.
 
For the differences between full MLE and lumping MLE, please refer to our article **Reverse Transcription Errors and RNA-DNA Differences at Short Tandem Repeats** by Arkarachai Fungtammasan, Marta Tomaszkiewicz, Rebeca Campos-Sanchez, Kristin Eckert, Michael DeGiorgio, and Kateryna D. Makova.
 
 
## System requirement
These scripts were written in R, a language and environment for statistical computing. The scripts also need `snow` package for multiprocessing and `combinat` package for combinatorial calculation.
 
 
## Installation
The script can be called directly given that the R program and all additional libraries were installed.
 
 
## Usage
To operate the script, run the command `Rscript  script_name  --filename=input_file_name`. The standard output will contain all the sets of parameter and their likelihood values that the algorithm were searched through. The last set of parameters will be the best estimated set of parameters given the initial set of parameters. Initial sets of parameters were chosen randomly in log scale.
To see the parameter option in help page, type `Rscript script_name` or `Rscript script_name --help`. This will print out all the parameter option that the user can change. Then the user can change the parameter through keyword argument. For example, to change the bin size in MLE calculation and number of core, the use can type:
`Rscript fullMLE_2lib.R --filename=inputfile.txt --bin_size=3 --N_core=4 > outputfile.txt`
 
 
## Data format
The input is in tab delimited format of seven columns Including: locus_name, genotype, STR_motif, STR_class, functional_compoment, STR_length_profile_RNAseq_lib1, STR_length_profile_ RNAseq_lib2. The script only uses of genotype, STR_motif, STR_length profile_RNAseq _lib1, and STR_length_profile_ RNAseq _lib2; therefore, the rest of value can be random character or dot. Each file should contain only STR that have the same motif and homozygote genotype. Both genotype and STR length profile of RNA-seq are coded as numbers of base pair. The STR lengths in STR length profiles are not need to be sorted.
 
Here is example of input
> chr14_901176_901185       	9 A 1 intron 10,9,9,9,9,9,9,9,9,9,9  	10,9,9,9,9,9,9,9

> chr14_901492_901501       	9 A 1 intron 10,9,9,9,9,9,9,9,9,9     	10,10,8,9,9,9,9,9,9,9


## Citing *STR-RNA-MLE*
Arkarachai Fungtammasan, Marta Tomaszkiewicz, Rebeca Campos-Sanchez, Kristin Eckert, Michael DeGiorgio, and Kateryna D. Makova,  Reverse Transcription Errors and RNA-DNA Differences at Short Tandem Repeats 
 
