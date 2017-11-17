##### README of MultivariateGWAMA.py #####
#
# 11 Oct 2017
# Kyoko Watanabe (k.watanabe@vu.nl)
##########################################

1. How to use
1.1 python dependency
Primary python version shoudl be 2.7
Make sure all the following packages are installed.
pandas, numpy, argparse, scipy and tempfile

1.2 Execution
The script can be executed by the following command
> python MultivariableGWAMA.py -c <config file> -i <intercept file>
Two flags, -c/--config and -i/--intercept are required. See below for format of the files.
Optionally -o/--out flag can be used to specify prefix of the output file. For example, "-o test" will create output file "test.txt". If the option is not given, default oputput file is "multivariateGWAS.txt".
--twoside is an option to align direction of alleles and convertion of Z <-> P are based on two sided. Otherwise, direction is not aligned with one side conversion.

1.3 input file format
GWAS summary statistics can be either plain text file or gzipped file with either tab or single space separated.
Config file is to specify column names, sample size and file name for each GWAS. See below for an example.
chrom, pos, a1, a2, weight/N and p columns are mandatory. a1 is the effect allele. When neither effect nor oddsratio is given, the program assumes a1 allele has increase effects on phenotype.

---------------------------------------
chrom chr # this means chromsome column name is "chr"
pos bp
a1 effect_allele
a2 non_effect_allele
effect beta #this should be signed effect size
weight N #when sample size is avilable in GWAS file
p p
rsID rsID
process GWAS1.txt

chrom chr
pos bp
a1 effect_allele
a2 non_effect_allele
oddratio or #if you only have or insterad effect size, you can use oddsratio flag
N 85000 #if per SNP sample size is not availabe, you can set total sample size here
p p
rsID rsID
process GWAS2.txt
---------------------------------------

Intercept is the gcov_int from LD score regression. To make sure, this is not actual sample overlap but takes phenotipic corration into account. The intercept C can be defined as Ns*rp/Na*Nb where Na and Nb are the sample size of GWAS a and b, Ns is the actual sample overlap, rp is the phenotypic correlation of GWAS a nad b.
The input file should contain lower triangle of pair-wise intercept matrix excluding diagonal. And should be ordered sample as the config file. For example, if you want to analyse GWAS a, b and c, the intercept file should look like

---------------------------------------
Cba
Cca Ccb
---------------------------------------

where Cij is the intercept between GWAS i and j, without header.

1.4 Output file
Column name should be explanatory except Neff. Neff is still under discussion so please do not use this value.
The weight column is the simple sum of the sample size so not considering sample overlap here.

2. Reference
Baselmans, B.M. et al. Multivariate genome-wide and integrated transcriptome and epigenome-wide analyses of the well-being spectrum. bioRxiv. https://dot.org/10.1101/115915
