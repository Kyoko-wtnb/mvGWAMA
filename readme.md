mvGWAMA: Multivariate GWAS meta-analysis
========================================
mvGWAMA is a python script to perform a GWAS meta-analysis when there are sample overlap.

current version: v.0.0.1  
last update: 2018-01-31

## Citation
In preparation

## Requirements
* ```Python 3> version >= 2.7```
* ```pandas version >= 0.19.0```
* ```numpy version >= 1.11.2```
* ```argparse version >= 1.2.1```
* ```scipy version >= 0.13.3```

## Getting started
You can either clone this repository or simply download python script.
```
git clone https://github.com/Kyoko-wtnb/mvGWAMA.git
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
* ```--twoside```: By default, direction of effect is not aligned and conversion between P-value and Z-score is one-sided.
When this flag is provided, direction is aligned and conversion between P-value and Z-score is two-sided.
This flag is highly recommended for a meta-analysis of the same or similar phenotypes.
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
The intercept (*CTI*) is defined as *N<sub>s</sub>* x *r<sub>p</sub>* / *N<sub>a</sub>N<sub>b</sub>* where *N<sub>a</sub>* and *N<sub>b</sub>* are the sample size of GWAS *a* and *b*, *N<sub>s</sub>* is the actual sample overlap, *r<sub>p</sub>* is the phenotypic correlation of phenotype *a* and *b*. However, sample overlap and phenotypic correlation are often not available. In that case, we can estimate intercept by using LD score regression [PMID:25642630] (the cross-trait LD score regression intercept).
The input file should contain lower triangle of pair-wise intercept matrix excluding diagonal. And the rows and columns should be ordered same as the order of GWAS file in the config file. For example, if you want to analyse GWAS *a*, *b* and *c*, the intercept file should look like

```
CTI_ba
CTI_ca CTI_cb
```
where *CTI<sub>ij</sub>* is the intercept between GWAS *i* and *j*, without header.

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
When this scripts is used on the HPC cluster with access to multiple cores, you can parallelize process of each chromosome as the following.
Note that, parallelizing 22 chromosome at once might cause memory error depends on how much memory access you have.
If memory error occur, please split chromosomes into multiple loops.

```
for i in {1..22}
do
(
	python mvGWAMA.py -c <config file> -i <intercept file> -ch $i -o temp$i.txt
)&
done
wait

for i in {2..22}
do
	awk 'NR>=2' temp$i.txt >>temp1.txt
done

mv temp1.txt <full output file>
```

## Effective sample size per SNP
The effective sample size (*N<sub>eff</sub>*) is computed for each SNP *k* from the matrix *M*, containing the sample size *N<sub>i</sub>* of each cohort *i* on the diagonal and the estimated number shared data points *N<sub>sij</sub>* x *rho<sub>ij</sub>* = *CTI<sub>ij</sub>* x sqrt(*N<sub>i</sub>N<sub>j</sub>*) for each pair of cohorts *i* and *j* as the off-diagonal values.
*N<sub>eff</sub>* is computed recursively as follows. Starting with the first cohort in *M*, *N<sub>eff</sub>* is first increased by *M<sub>1,1</sub>*, corresponding to the sample size of that cohort.
The proportion of samples shared between cohort 1 and each other cohort *j* is them computed as *p<sub>1,j</sub>* = *M<sub>1,j</sub>* / *M<sub>j,j</sub>*, and *M* is then adjusted to remove this overlap, multiplying all values in each column *j* by 1-*p<sub>1,j</sub>*.
This amounts to reducing the sample size of each other cohort *j* by the number of samples it shares with cohort 1, and reducing the shared samples between cohort *j* and subsequent cohorts by the same proportion.
After this, the first row and column of *M* are discarded, and the same process is applied to the new *M* matrix. This is repeated until *M* is empty.

Note that computing per SNP effect size takes longer than without --neff-per-snp flag.

```
function compute_neff(M):
	n_eff = M[1,1]
	if (dimension(M) > 1) n_eff += compute_neff(reduce(M))
	return n_eff

function reduce(M):
	for (j := 2 to dimension(M))
		p = M[1,j]/M[j,j]
		for (i := 2 to dimension(M))
			M[i,j] = (1-p)*M[i,j]
	return M[-1,-1]
```

## Note and tips
### 1. rsID
When rsID column is provided for any of input GWAS files, rsID is extracted from those files.
For SNPs whoes rsID is not available, the rsID column is substituted with unique ID of the SNP which consists of chr:position:alleleA_alleleB where alleles are alphabetically ordered.
For SNPs that have rsID in more than one input file but rsID does not match across files, the rsID column is substituted with unique ID.
Note that this script does not check the duplication of rsID.

### 2. Alignment of direction
Effect allele is defined as the one provided in the first GWAS file.
Therefore, the alleles of the 2nd or later GWAS files are aligned to the ones in the first file.

### 3. Duplicated SNPs
Prior to the analysis, duplicated SNPs based on the unique ID (chr:pos:alleleA_alleleB where alleles are alphabetically ordered) are removed.

### 4. Matching SNPs between GWAS
SNPs are matched based on unique ID (chr:pos:alleleA_alleleB where alleles are alphabetically ordered) between input GWAS files.
Therefore, please make sure that all of input GWAS files are based on the same genome reference (e.g. hg19).

## Updates
2018-01-31: v0.0.1 (First release)

## Licence
This script is licensed under GNU GPL v3.

## Authors
Kyoko Watanabe (VU University Amsterdam)  
Christiaan de Leeuw (VU University Amsterdam)
