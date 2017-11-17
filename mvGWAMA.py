#!/usr/bin/python

'''
Multivariate genome-wide assocaition meta analysis

8 Nov 2017
Kyoko Watanabe (k.watanabe@vu.nl)
'''

import sys
import os
import pandas as pd
import numpy as np
import re
import argparse
import scipy.stats as st
import math
import time
from tempfile import mkdtemp
from joblib import Parallel, delayed

__version__ = '0.0.0'
__date__ = '17/Nov/2017'
HEADMSS = "#####################################################\n"
HEADMSS += "# Multivariate genome-wide association meta-analysis\n"
HEADMSS += "# Version: {V}\n".format(V=__version__)
HEADMSS += "# Last update: {D}\n".format(D=__date__)
HEADMSS += "# (c) 2017 Kyoko Watanabe\n"
HEADMSS += "# GNU General Public Licence v3\n"
HEADMSS += "#####################################################\n"

parser = argparse.ArgumentParser()
parser.add_argument('-c', '--config', default=None, type=str, help="(Required) Config file of summary statistics.")
parser.add_argument('-i', '--intercept', default=None, type=str, help="(Required) File name of the intercelpt matrix (lower triangle).")
parser.add_argument('-o', '--out', default="mvGWAMA", type=str, help="Output file name. 'mvGWAMA' by default.")
parser.add_argument('-ch', '--chrom', default=None, type=int, help="To run for a specific chromosome.")
parser.add_argument('-p', '--parallel', default=None, type=int, help="To parallelize process, provide the number of cores/thread.")
parser.add_argument('--twoside', default=False, action='store_true', help="Use this flag to convert P to Z by two sided.")
parser.add_argument('--neff-per-snp', default=False, action='store_true', help="Use this flag to compute effective samplesize per SNP (runtime will be longer). Otherwise, per SNP effect size is computed based on proportion of total Neff to total Nsum.")
parser.add_argument('--no-weight', default=False, action='store_true', help="Use this flag to not weight by sample size.")

### global variable
tmpdir = os.path.join(mkdtemp()) #path to temp for memmap files
allele_idx = ['A', 'C', 'G', 'T']
allele_map = {'A':0, 'C':1, 'G':2, 'T':3}
nGWAS = 0
nSNPs = [0]*23
Nall = []

##### Return index of a1 which exist in a2 #####
def ArrayIn(a1, a2):
	results = np.where(np.in1d(a1, a2))[0]
	return results

##### Return index of a1 which do not exist in a2 #####
def ArrayNotIn(a1, a2):
    tmp = np.where(np.in1d(a1, a2))[0]
    return list(set(range(0,len(a1)))-set(tmp))

##### return unique element in a list #####
def unique(a):
	unique = []
	[unique.append(s) for s in a if s not in unique]
	return unique

##### return index of duplicated values #####
def duplicated(a):
	unq, unq_count = np.unique(a, return_counts=True)
	duplicated_value = unq[unq_count > 1]
	duplicated_idx = ArrayIn(a, duplicated_value)
	return duplicated_idx

##### return index of non-duplicated values #####
def non_duplicated(a):
	unq, unq_count = np.unique(a, return_counts=True)
	return ArrayIn(a, unq[unq_count == 1])

### read intercept matrix
# input file should contain lower triangle (excluding diagonal)
# take absolute value (v 0.0.0)
# C is winsolized between -1 and 1
def getIntercept(infile):
	C = []
	with open(infile, 'r') as fin:
		for l in fin:
			l = l.strip().split()
			tmp = [abs(float(x)) if float(x)>=-1 and float(x)<=1 else float(1) for x in l]
			C.append(tmp)
	return C

def countGWASfiles(infile):
	n = 0
	with open(infile, 'r') as inf:
		for l in inf:
			if l.startswith("process"):
				n += 1
				f = l.split()[1]
				if not os.path.isfile(f):
					sys.exit("\n ERROR: Input GWAS file '"+f+"' does not exist.")
	return n

### match rsID
def match_rsID(ids):
	if ids[0] == "NA":
		return ids[1]
	elif ids[1] == "NA":
		return ids[0]
	elif ids[0] == ids[1]:
		return ids[0]
	else:
		return "NA"

### Update matrices
def updateMatrix(gwas, chrom, GWASidx, C, nsnps, noweight):
	# global nSNPs
	if GWASidx == 1 or not os.path.isfile(tmpdir+'/snps'+str(chrom)+'.dat'):
		### initialize snps, info, w and v
		# snps: 0:chr, 1:pos, 2:a1, 3:a2 (int)
		# info: 0:uniqID, 1:rsID, 2:direction (str)
		# w : store weight per GWAS (int)
		# v: 0:wz, 1:w2, 2:wwc

		print "Initializing snps matrix..."
		snps = np.memmap(tmpdir+'/snps'+str(chrom)+'.dat', dtype='int64', mode='w+', shape=(len(gwas), 4))
		snps[:] = gwas[:,0:4]
		snps.flush()
		info = np.memmap(tmpdir+'/info'+str(chrom)+'.dat', dtype='S'+str(nGWAS+10), mode='w+', shape=(len(gwas), 3))
		info[:] = np.c_[[str(l[1])+":"+"_".join(sorted([str(l[2]), str(l[3])])) for l in gwas], gwas[:,7], ["?"*(GWASidx-1)+'+' if x>0 else "?"*(GWASidx-1)+'-' for x in gwas[:,5].astype(int)]]
		info.flush()

		print "Initializing weight matrix..."
		w = np.memmap(tmpdir+'/w'+str(chrom)+'.dat', dtype='int64', mode='w+', shape=(len(gwas), nGWAS))
		w[:,0] = gwas[:,6]
		w.flush()

		print "Initializing variable matrix..."
		v = np.memmap(tmpdir+'/v'+str(chrom)+'.dat', dtype='float128', mode='w+', shape=(len(gwas), 3))
		if noweight:
			v[:,0] = gwas[:,4]
			v[:,1] = 1
		else:
			v[:,0] = np.multiply(np.sqrt(gwas[:,6].astype(int)), gwas[:,4])
			v[:,1] = gwas[:,6]
		v[:,2] = 0
		v.flush()

		nsnps = len(snps)
		del snps, info, w, v
	else:
		print "Checking additional SNPs..."
		info = np.memmap(tmpdir+'/info'+str(chrom)+'.dat', dtype='S'+str(nGWAS+10), mode='r', shape=(nsnps, 3), order='C')
		cur_uid = [str(l[1])+":"+"_".join(sorted([str(l[2]), str(l[3])])) for l in gwas]
		new_idx = ArrayNotIn(cur_uid, info[:,0])
		print "Detected "+str(len(new_idx))+" additional SNPs"
		del info

		print "Loading matrices..."
		snps = np.memmap(tmpdir+'/snps'+str(chrom)+'.dat', dtype='int64', mode='r+', shape=(nsnps+len(new_idx), 4), order='C')
		info = np.memmap(tmpdir+'/info'+str(chrom)+'.dat', dtype='S'+str(nGWAS+10), mode='r+', shape=(nsnps+len(new_idx), 3), order='C')
		w = np.memmap(tmpdir+'/w'+str(chrom)+'.dat', dtype='int64', mode='r+', shape=(nsnps+len(new_idx), nGWAS), order='C')
		v = np.memmap(tmpdir+'/v'+str(chrom)+'.dat', dtype='float128', mode='r+', shape=(nsnps+len(new_idx), 3), order='C')

		n = ArrayIn(cur_uid, info[:,0])
		m = ArrayIn(info[:,0], cur_uid)

		print "Aligning direction..."
		gwas[n,5] = [l[2] if l[0]==l[1] else -1*l[2] for l in np.c_[snps[m,2], gwas[n,2], gwas[n,5]]]

		print "Updating matrices..."
		snps[nsnps:] = gwas[new_idx, 0:4]
		info[m,1] = [match_rsID(l) for l in np.c_[info[m,1], gwas[n,7]]]
		info[m,2] = [l[1]+'+' if l[0]>0 else l[1]+'-' for l in np.c_[gwas[n,5], info[m,2]]]
		info[nsnps:] = np.c_[[str(l[1])+":"+"_".join(sorted([str(l[2]), str(l[3])])) for l in gwas[new_idx]],
			gwas[new_idx,7],
			["?"*(GWASidx-1)+'+' if x>0 else "?"*(GWASidx-1)+'-' for x in gwas[new_idx,5].astype(int)]]
		info[ArrayNotIn(info[0:nsnps,0], cur_uid),2] = [x+"?" for x in info[ArrayNotIn(info[0:nsnps,0], cur_uid),2]]
		w[m,GWASidx-1] = gwas[n,6]
		w[nsnps:,GWASidx-1] = gwas[new_idx,6]

		if noweight:
			v[m,0] = np.add(v[m,0], gwas[n,4])
			v[m,1] = np.add(v[m,1], [1]*len(m))
			for i in range(1,GWASidx):
				v[m,2] = np.add(v[m,2], C[GWASidx-2][i-1])
			v[nsnps:,0] = gwas[new_idx,4]
			v[nsnps:,1] = 1
			v[nsnps:,2] = 0
		else:
			v[m,0] = np.add(v[m,0], np.multiply(np.sqrt(gwas[n,6].astype(int)), gwas[n,4]))
			v[m,1] = np.add(v[m,1], gwas[n,6])
			for i in range(1,GWASidx):
				v[m,2] = np.add(v[m,2], np.multiply(np.sqrt(w[m,i-1]),np.sqrt(w[m,GWASidx-1]))*C[GWASidx-2][i-1])
			v[nsnps:,0] = np.multiply(np.sqrt(gwas[new_idx,6].astype(int)), gwas[new_idx,4])
			v[nsnps:,1] = gwas[new_idx,6]
			v[nsnps:,2] = 0

		### sort by position
		n = snps[:,1].argsort()
		snps[:] = snps[n,:]
		info[:] = info[n,:]
		w[:] = w[n,:]
		v[:] = v[n,:]
		nsnps = len(snps)

		del snps, info, w, v
	return nsnps

### Process each GWAS sumstat file
def processFile(gwasfile, C, GWASidx, chrom, pos, a1, a2, p, effect, oddsratio, N, weight, rsID, delim, args):
	global allele_idx
	global allele_map
	global nSNPs
	global Nall

	cols = [chrom, pos, a1, a2, p]
	if rsID is not None:
		cols.append(rsID)
	if effect is not None:
		cols.append(effect)
	if oddsratio is not None:
		cols.append(oddsratio)
	if weight is not None:
		cols.append(weight)
	gwas = pd.read_table(gwasfile, sep=delim, header=0, usecols=cols)
	header = list(gwas.columns.values)
	gwas = np.array(gwas)
	if type(gwas[0,header.index(chrom)]) is str:
		gwas[:,header.index(chrom)] = [x.replace('chr','').replace('X','23').replace('x','23') for x in gwas[:,header.index(chrom)]]

	if N is not None:
		Nall.append(int(N))
	else:
		Nall.append(max(gwas[:,header.index(weight)].astype(int)))

	### filter on chr
	if args.chrom is not None:
		print "Filtering on chromosome "+str(args.chrom)
		gwas = gwas[gwas[:,header.index(chrom)].astype(int)==args.chrom]

	print "Detected "+str(len(gwas))+" SNPs in the file"
	### check header

	### check effect
	# if oddsratio is given instead effec, take log
	# if both are not given, assume a1 alele has increasing risk
	# convert effect to 1/-1
	print "Checking effect column..."
	if effect is not None:
		effect = [1 if x>0 else -1 for x in gwas[:,header.index(effect)]]
	elif oddsratio is not None:
		effect = [1 if x>1 else -1 for x in gwas[:,header.index(oddsratio)]]
	else:
		print "WARNING: Neither signed effect size or odds ration was gievn, a1 allele is considered as risk increasing allele."
		effect = [1 for i in range(0, len(gwas))]

	### check weight
	# if weight is not given, assign N to all SNPs
	print "Checking weight column..."
	if weight is not None:
		weight = gwas[:,header.index(weight)].astype(int)
	else:
		weight = [int(N)]*len(gwas)

	### check rsID
	# if rsID is not given, store "NA"
	print "Checking rsID column..."
	if rsID is not None:
		rsID = gwas[:,header.index(rsID)]
	else:
		rsID = ["NA"]*len(gwas)

	### reformat gwas
	# 0:chr, 1:pos, 2:a1, 3:a2, 4:p(Z later), 5:effect, 6:weight, 7:rsID
	print "Formatting gwas input..."
	gwas = gwas[:, [header.index(chrom), header.index(pos), header.index(a1), header.index(a2), header.index(p)]]
	# allele convert to int
	tmp_a = unique(gwas[:,2])
	tmp_a = [x.upper() for x in unique(tmp_a+gwas[:,3].tolist())]
	for a in tmp_a:
		if a not in allele_map:
			allele_idx.append(a)
			allele_map[a] = len(allele_idx)-1
	gwas = np.c_[gwas[:,[0,1]], [allele_map[x.upper()] for x in gwas[:,2]], [allele_map[x.upper()] for x in gwas[:,3]], gwas[:,4:], effect, weight, rsID]
	effect = None
	weight = None
	rsID = None
	tmp_a = None

	### remove duplicated SNPs
	n = non_duplicated([str(l[0])+":"+str(l[1])+":"+"_".join(sorted([str(l[2]), str(l[3])])) for l in gwas[:,0:4]])
	if len(n) < len(gwas):
		print "Warning: "+str(len(gwas)-len(n))+" SNPs are removed due to duplicated uniqID."
		gwas = gwas.take(n,0)

	### sort gwas by chr and pos
	n = np.lexsort((gwas[:,0].astype(int), gwas[:,1].astype(int)))
	gwas = gwas.take(n,0)

	### compute Z
	print "Converting P to Z score..."
	# replace P == 0 to the minimum P-value in the input file
	if len(np.where(gwas[:,4]==0.0)[0])>0:
		print "WARNING: P-value < 1e-323 is replaced with 1e-323"
		gwas[gwas[:,4]==0.0,4] = 1e-323
	if len(np.where(gwas[:,4]==1)[0])>0:
		print "WARNING: P-value 1 is replaced with 0.999999"
		gwas[gwas[:,4]==1,4] = 0.999999
	if args.twoside:
		gwas[:,4] = -1.0*gwas[:,5]*st.norm.ppf(list(np.divide(gwas[:,4],2)))
	else:
		gwas[:,4] = -1.0*st.norm.ppf(list(gwas[:,4]))

	chroms = unique(gwas[:,0])
	## parallelize
	if args.parallel is not None and args.parallel>0:
		if args.parallel > len(chroms):
			args.parallel = len(chroms)
		# tmp_nSNPs = Parallel(n_jobs=args.parallel)(delayed(updateMatrix)(gwas[gwas[:,0]==c], c, GWASidx, C, nSNPs[c-1], args.no_weight) for c in chroms)
		for i in range(0,len(chroms)):
			nSNPs[chroms[i]-1] = tmp_nSNPs[i]
	else:
		for c in chroms:
			nSNPs[int(c)-1] = updateMatrix(gwas[gwas[:,0]==c], int(c), GWASidx, C, nSNPs[int(c)-1], args.no_weight)


### process GWAS
def processGWAS(C, args):
	GWASidx = 0
	snps = None
	chrom = None
	pos = None
	a1 = None
	a2 = None
	rsID = None
	p = None
	effect = None
	oddsratio = None
	N = None
	weight = None
	delim = "\t"

	with open(args.config, 'r') as inconfig:
		for l in inconfig:
			if l=="\n":
				continue
			if re.match(r'^#', l):
				continue
			l = l.strip().split()
			if l[0] == "chrom":
				chrom = l[1]
			elif l[0] == "pos":
				pos = l[1]
			elif l[0] == "a1":
				a1 = l[1]
			elif l[0] == "a2":
				a2 = l[1]
			elif l[0] == "rsID":
				rsID = l[1]
			elif l[0] == "p":
				p = l[1]
			elif l[0] == "effect":
				effect = l[1]
			elif l[0] == "oddsratio":
				oddsratio = l[1]
			elif l[0] == "N":
				N = int(l[1])
			elif l[0] == "weight":
				weight = l[1]
			elif l[0] == "delim":
				delim = l[1]
			elif l[0] == "process":
				gwasfile = l[1]
				GWASidx += 1
				print "------------------------------------------------"
				print "Process GWAS "+str(GWASidx)+": "+gwasfile
				if not (chrom and pos and a1 and a2 and p):
					sys.exit("\nERROR: Not enought columns are provided in the config file. Chrom, pos, a1, a2 and p columns are required for all input GWAS files.")
				if not (N or weight):
					sys.exit("\nERROR: Neither N nor weight are provided in the config file.")
				processFile(gwasfile, C, GWASidx, chrom, pos, a1, a2, p, effect, oddsratio, N, weight, rsID, delim, args)
				chrom = None
				pos = None
				a1 = None
				a2 = None
				rsID = None
				p = None
				effect = None
				oddsratio = None
				N = None
				weight = None

		# v = np.c_[v, [int(round(sum(np.square(x)))) for x in w[:,1:]]]
		return

##### compute z from stored variables and combert to P #####
def computeZ(twoside):
	out = []
	for chrom in range(1,24):
		if os.path.isfile(tmpdir+'/snps'+str(chrom)+'.dat'):
			snps = np.memmap(tmpdir+'/snps'+str(chrom)+'.dat', dtype='int64', mode='r+', shape=(nSNPs[chrom-1], 4), order='C')
			info = np.memmap(tmpdir+'/info'+str(chrom)+'.dat', dtype='S'+str(nGWAS+10), mode='r+', shape=(nSNPs[chrom-1], 3), order='C')
			w = np.memmap(tmpdir+'/w'+str(chrom)+'.dat', dtype='int64', mode='r+', shape=(nSNPs[chrom-1], nGWAS), order='C')
			v = np.memmap(tmpdir+'/v'+str(chrom)+'.dat', dtype='float128', mode='r+', shape=(nSNPs[chrom-1], 3), order='C')
			z = np.divide(v[:,0].astype(float), np.sqrt(np.add(v[:,1].astype(float), 2*v[:,2].astype(float))))
			if twoside:
				p = st.norm.cdf(list(-1.0*np.absolute(z)))*2
			else:
				p = st.norm.cdf(list(-1.0*z))
			if len(out)==0:
				out = np.c_[snps[:,[0,1]], [allele_idx[x] for x in snps[:,2]], [allele_idx[x] for x in snps[:,3]], info[:,1], z, p, [sum(x) for x in w], info[:,2]]
			else:
				out = np.r_[out, np.c_[snps[:,[0,1]], [allele_idx[x] for x in snps[:,2]], [allele_idx[x] for x in snps[:,3]], info[:,1], z, p, [sum(x) for x in w], info[:,2]]]
	### return chr, pos, a1, a2, rsID, z, p, weight, direction
	return out

##### reduce N matrix
def reduceMat(M):
	if M[0,0]<0: return M[1:,1:]
	prop = M[1:,0]/np.diag(M)[1:]
	return M[1:,1:]*(1-prop)

##### compute effective N recursively
def NeffMap(M):
	if M[0,0]<0: M[0,0]=0
	if len(M)<=1: return M[0,0]
	else: return M[0,0]+NeffMap(reduceMat(M))

def getNeffPerSNP(C):
	Neff = []
	for chrom in range(1,24):
		if os.path.isfile(tmpdir+'/snps'+str(chrom)+'.dat'):
			w = np.memmap(tmpdir+'/w'+str(chrom)+'.dat', dtype='int64', mode='r+', shape=(nSNPs[chrom-1], nGWAS), order='C')
			for l in w:
				Nmat = []
				for i in range(0,nGWAS):
					if i==0:
						Nmat.append([l[i]]+[0]*(nGWAS-1))
					else:
						tmp = []
						for j in range(0,i):
							tmp.append(math.sqrt(l[i]*l[j])*C[i-1][j])
						Nmat.append(tmp+[l[i]]+[0]*(nGWAS-1-i))
				Nmat = np.array(Nmat).astype(float)
				n = np.where(np.diag(Nmat)>0)
				Neff.append(NeffMap(Nmat[n,n]))
	return Neff

def getNeff(C):
	global Nall
	Nmat = []
	for i in range(0,nGWAS):
		if i==0:
			Nmat.append([Nall[i]]+[0]*(nGWAS-1))
		else:
			tmp = []
			for j in range(0,i):
				tmp.append(math.sqrt(Nall[i]*Nall[j])*C[i-1][j])
			Nmat.append(tmp+[Nall[i]]+[0]*(nGWAS-1-i))
	Nmat = np.array(Nmat).astype(float)
	return NeffMap(Nmat)

def main(args):
	start_time = time.time()

	global HEADMSS
	HEADMSS += "Flags used:\n"
	opts = vars(args)
	options = ['--'+x.replace('_', '-')+' '+str(opts[x])+' \\' for x in opts.keys() if opts[x]]
	HEADMSS += "\t"+"\n\t".join(options).replace('True','').replace('False','')
	print HEADMSS

	### check arguments
	if args.config is None:
		parser.print_help()
		sys.exit("\nERROR: Config file is required.")
	if args.intercept is None:
		parser.print_help()
		sys.exit("\nERROR: Intercept file is required.")

	### check input files
	if not os.path.isfile(args.config):
		sys.exit("\nERROR: Config file '"+args.config+"' does not exist.")
	if not os.path.isfile(args.intercept):
		sys.exit("\nERROR: Intercept file '"+args.intercept+"' does not exist.")

	### count the number of GWAS to process and check if the file exist
	global nGWAS
	nGWAS = countGWASfiles(args.config)
	print "\nDetected "+str(nGWAS)+" input GWAS files.\n"

	### get intercept matrix
	C = getIntercept(args.intercept)
	if len(C) != nGWAS-1:
		sys.exit("\nERROR: The dimention of intercept matrix is wrong. The matrix should be lower off diagonal of pariwise intercept.")

	### process files and store variables
	processGWAS(C, args)
	print "------------------------------------------------\n"

	### compute test statistics
	results = computeZ(args.twoside)
	## replace rsID=NA to uniqID
	results[results[:,4]=="NA",4] = [[str(l[0])+":"+str(l[1])+":"+"_".join(sorted([l[2], l[3]])) for l in results[results[:,4]=="NA",0:4]]]

	### compute effective sample size
	Neff_total = getNeff(C)
	Nprop = Neff_total/sum(Nall)
	print "***Sample size***"
	print "Sum of sample size: "+str(sum(Nall))
	print "Effective sample size: "+str(round(Neff_total,2))
	print "Proportion of Neff to Nsum: "+str(round(Nprop,4))

	if args.neff_per_snp:
		print "Computing per SNP Neff..."
		Neff = getNeffPerSNP(C)
		results = np.c_[results[:,0:8], [round(x,2) for x in Neff], results[:,8]]
	else:
		print "WARNING: Use ratio of total Neff to total Nsum to compute per SNP Neff"
		print "         Use --neff-per-snp flag to compute accurate per SNP Neff."
		results = np.c_[results[:,0:8], [round(x,2) for x in results[:,7].astype(float)*Nprop], results[:,8]]

	with open(args.out+".txt", 'w') as o:
		o.write("\t".join(["chr", "pos", "a1", "a2", "rsID", "z", "p", "Nsum", "Neff", "dir"])+"\n")
	with open(args.out+".txt", 'a') as o:
		np.savetxt(o, results, delimiter="\t", fmt="%s")
	os.system("rm -r "+tmpdir)
	print "\nProcess completed\nProgram run time: "+str(round(time.time()-start_time,2))+" sec"

if __name__ == "__main__": main(parser.parse_args())
