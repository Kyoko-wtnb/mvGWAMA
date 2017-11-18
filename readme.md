mvGWAMA: Multivariate GWAS meta-analysis
========================================
mvGWAMA is a python script to perform a GWAS meta-analysis when there are sample overlap.

current version: v.0.0.0  
last update: 17/Nov/2017

## Citation
under preparation

## Requirements
* ```Python 3> version >= 2.7```
* ```pandas version >= 0.19.0```
* ```numpy version >= 1.11.2```
* ```argparse version >= 1.2.1```
* ```scipy version >= 0.13.3```

## Getting started
You can either clone this repository or simply download python script.
```
git clone
```

The script can be executed by the following command.
There are two required arguments.
```
python mvGWAMA.py -c <config file> -i <intercept file>
```

In addition, there are several optional arguments.
```
optional arguments:
  -h, --help            show this help message and exit
  -c CONFIG, --config CONFIG
                        (Required) Config file of summary statistics.
  -i INTERCEPT, --intercept INTERCEPT
                        (Required) File name of the intercelpt matrix (lower
                        triangle).
  -o OUT, --out OUT     Output file name. 'mvGWAMA' by default.
  -ch CHROM, --chrom CHROM
                        To run for a specific chromosome.
  -p PARALLEL, --parallel PARALLEL
                        To parallelize process, provide the number of
                        cores/thread.
  --twoside             Use this flag to convert P to Z by two sided with
  						alignment of direction of effects.
  --neff-per-snp        Use this flag to compute effective samplesize per SNP
                        (runtime will be longer). Otherwise, per SNP effect
                        size is computed based on proportion of total Neff to
                        total Nsum.
  --no-weight           Use this flag to not weight by sample size.
```

* ```-c / --config```: (Required) Path to a input config file (see below for format).
* ```-i / --intercept```: (Required) Path to a input intercept file (see below for format).
* ```-o / --out```: Prefix of the output file. Default is mvGWAMA which creates two output files, mvGWAMA.txt and mvGWAMA.log.
* ```-ch / --chrom```: To perform mvGWAMA for a specific chromosome, this flag can be used with integer variable between 1 and 23.
* ```-p / --parallel```: To parallelize a process, a number of core can be provided with this flag.
See section xxx for details of parallelization.
* ```--twoside```: By default, direction of effect is not aligned and conversion between P-value and Z-score is one-sided.
When this flag is provided, direction is aligned and conversion betweeen P-value and Z-score is two-sided.
This flag is highhly recommended for a meta-analysis of the same or similar phenotypes.
* ```--neff-per-snp```: When this flag is provided, effective sample size is computed per SNP.
Otherwise, per SNP effective sample size is based on the ration of total effective sample size to total sum of sample size.
The total run time can be longer (around 20%) when this flag is used. See section xxx for details.
* ```--no-weight```: When this flag is provided, the weight is replaced by 1, otherwise, per SNP sample size is used as weight.

## Input file format
### 1. GWAS summary statistics and config file
GWAS summary statistics can be either plain text file or gzipped file with either tab or single space separated.
Config file is to specify column names, sample size and file name for each GWAS. See below for an example.  
chrom, pos, a1, a2, weight/N and p columns are mandatory. a1 is the effect allele. When neither effect nor oddsratio is given, the program assumes a1 allele has increase effects on phenotype.
rsID is not mandatory, and if it is not provided in any of GWAS files, unique ID is provided in the output file.  
Each GWAS file needs to be specified with "process" flag and all the column names have to be specified before that column.  
The lines start with "#" will be ignored.

```
chrom chr # this means chromsome column name is "chr"
pos bp
a1 effect_allele
a2 non_effect_allele
effect beta #this should be signed effect size
weight N #when sample size is available in GWAS file
p p
rsID rsID
process GWAS1.txt

chrom chr
pos bp
a1 effect_allele
a2 non_effect_allele
oddratio or #if you only have OR instead signed effect size, you can use oddsratio flag
N 85000 #if per SNP sample size is not available, you can set total sample size with N flag
p p
rsID rsID
process GWAS2.txt
```

### 2. intercept matrix
The intercept is defined as Ns*rp/Na*Nb where Na and Nb are the sample size of GWAS a and b, Ns is the actual sample overlap, rp is the phenotypic correlation of phenotype a and b. However, sample overlap and phenotypic correlation are often not available. In that case, we can estimate intercept by using LD score regression [ref] (the cross-trait LD score regression intercept).
The input file should contain lower triangle of pair-wise intercept matrix excluding diagonal. And the rows and columns should be ordered same as the order of GWAS file in the config file. For example, if you want to analyse GWAS a, b and c, the intercept file should look like

```
Cba
Cca Ccb
```
where Cij is the intercept between GWAS i and j, without header.

## Output file format
### 1. test statistics
The output text file contains the following 10 columns.

chr: chromsome  
pos: position  
a1: effect allele  
a2: other allele  
rsID: rsID  
z: z statistics  
p: P-value  
Nsum: sum of sample size of the SNP  
Neff: effective sample size of the SNP  
dir: directions of the SNP with the same order as the input GWAS files

### 2. log file
All the messages displayed on standard output is stored in the log file.
The total sum of sample size and total effective sample size are available in the bottom of the log file.

```
***Sample size***
Sum of sample size: 390469
Effective sample size: 384293.15
Proportion of Neff to Nsum: 0.9842
```

## Parallelization

## Effective sample size

## Note and tips
1. rsID
2. Alignment of direction
3. P-value
4. Duplicated SNPs
5. Matching SNPs between GWAS
6. Intercept
7. Effective sample size

## Updates

## Licence
This project is licensed under GNU GPL v3.

## Authors
Kyoko Watanabe (VU University Amsterdam)  
Christiaan de Leeuw (VU University Amsterdam)
