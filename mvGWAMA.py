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

__version__ = '0.0.0'
HEADMSS = "#####################################################\n"
HEADMSS += "# Multivariate genome-wide association meta-analysis\n"
HEADMSS += "# mvGWAMA.py\n"
HEADMSS += "# Version {V}\n".format(V=__version__)
HEADMSS += "# (c) 2017 Kyoko Watanabe\n"
HEADMSS += "# MIT Licence\n"
HEADMSS += "#####################################################\n"

parser = argparse.ArgumentParser()
parser.add_argument('-c', '--config', default=None, type=str, help="(Required) Config file of summary statistics.")
parser.add_argument('-i', '--intercept', default=None, type=str, help="(Required) File name of the intercelpt matrix (lower triangle).")
parser.add_argument('-o', '--out', default="multivariateGWAS", type=str, help="Output file name. 'multivariateGWAS' by default.")
parser.add_argument('--twoside', default=False, action='store_true', help="Use this flag to convert P to Z by two sided.")
parser.add_argument('-ch', '--chrom', default=None, type=int, help="To run for a specific chromosome.")
#parser.add_argument('--no-weight', default=False, action='store_true', help="Use this flag to not weight by sample size.")

### global variable
tmpdir = os.path.join(mkdtemp()) #path to temp for memmap files
allele_idx = ['A', 'C', 'G', 'T']
allele_map = {'A':0, 'C':1, 'G':2, 'T':3}

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
# C = [-1,1]
def getIntercept(infile):
	C = []
	with open(infile, 'r') as inf:
		for l in inf:
			l = l.strip().split()
			tmp = [abs(float(x)) if float(x)>=-1 and float(x)<=1 else float(1) for x in l]
			C.append(tmp)
	return C

### match rsID
def match_rsID(ids):
	if ids[1] == "NA":
		return ids[2]
	elif ids[2] == "NA":
		return ids[1]
	elif ids[1] == ids[2]:
		return ids[1]
	else:
		return ids[0]

### procedd each GWAS sumstat file
def processFile(gwasfile, C, snps, GWASidx, chrom, pos, a1, a2, p, effect, oddsratio, N, weight, rsID, delim, twoside, chromfilt):
	print "Process GWAS "+str(GWASidx)+": "+gwasfile
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

	### filter on chr
	if chromfilt is not None:
		print "Filtering on chromosome "+str(chromfilt)
		gwas = gwas[gwas[:,header.index(chrom)].astype(int)==chromfilt]

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
		effect = [1 for i in range(0, len(gwas))]

	### check weight
	# if weight is not given, assign N to all SNPs
	print "Checking weight column..."
	if weight is not None:
		weight = np.sqrt(gwas[:,header.index(weight)].astype(int))
	else:
		weight = np.full((len(gwas)), math.sqrt(N))

	### check rsID
	# if rsID is not given, store "NA"
	print "Checking rsID column..."
	if rsID is not None:
		rsID = gwas[:,header.index(rsID)]
	else:
		rsID = ["NA"]*len(gwas)

	### reformat gwas
	# 0:uniqID, 1:chr, 2:pos, 3:a1, 4:a2, 5:p(Z later), 6:effect, 7:weight, 8:rsID
	print "Formatting gwas input..."
	gwas = gwas[:, [header.index(chrom), header.index(pos), header.index(a1), header.index(a2), header.index(p)]]
	gwas[:,2] = np.char.upper(gwas[:,2].tolist())
	gwas[:,3] = np.char.upper(gwas[:,3].tolist())
	gwas = np.c_[[str(l[0])+":"+str(l[1])+":"+"_".join(sorted([l[2], l[3]])) for l in gwas[:,0:4]], gwas, effect, weight, rsID]
	effect = None
	weight = None
	rsID = None

	### remove duplicated SNPs
	n = non_duplicated(gwas[:,0])
	if len(n) < len(gwas):
		print "Warning: "+str(len(gwas)-len(n))+" SNPs are removed due to duplicated uniqID."
		#gwas = np.take(gwas, n)
		gwas = gwas.take(n,0)
	### sort gwas by uniqID
	n = gwas[:,0].argsort()
	gwas = gwas.take(n,0)

	### align allele and direction
	if GWASidx > 1:
		print "Aligning Allele"
		n = ArrayIn(gwas[:,0],snps[:,0])
		gwas[n,6] = [x[2] if x[0]==x[1] else -1*x[2] for x in np.c_[gwas[n,3], snps[ArrayIn(snps[:,0],gwas[:,0]),1], gwas[n,6]]]
	### compute Z
	print "Converting P to Z score..."
	# replace P == 0 to the minimum P-value in the input file
	gwas[gwas[:,5]==0.0,5] = 1e-324
	gwas[gwas[:,5]==1,5] = 0.999999
	if twoside:
		gwas[:,5] = -1.0*gwas[:,6]*st.norm.ppf(list(np.divide(gwas[:,5],2)))
	else:
		gwas[:,5] = -1.0*st.norm.ppf(list(gwas[:,5]))

	if GWASidx == 1:
		### initialize snps, w and v
		# snps: store, uniqID, a1 and a2
		# w : store weight per GWAS (squared N)
		# v: stroe each variables
		#	0:wz, 1:w2, 2:wwc, 3:direction, 4: rsID
		print "Initializing snps matrix..."
		snps = np.c_[gwas[:,[0,3,4]], ["+" if x>0 else "-" for x in gwas[:,6]], gwas[:,8]]
		print "Initializing weight matrix..."
		w = np.memmap(tmpdir+'/w.dat', dtype='float128', mode='w+', shape=(len(gwas), 1))
		w[:,0] = gwas[:,7]
		w.flush()
		print "Initializing variable matrix..."
		v = np.memmap(tmpdir+'/v.dat', dtype='float128', mode='w+', shape=(len(gwas), 3))
		v[:,0] = np.multiply(gwas[:,7], gwas[:,5])
		v[:,1] = np.square(gwas[:,7])
		v[:,2] = [0]*len(gwas)
		# v = np.c_[np.multiply(gwas[:,7], gwas[:,5]), np.square(gwas[:,7]), [0]*len(gwas)]
		v.flush()
	else:
		print "Checking additional SNPs..."
		new_uid = gwas[ArrayNotIn(gwas[:,0],snps[:,0]),0]
		w_old = np.memmap(tmpdir+'/w.dat', dtype='float128', mode='r', shape=(len(snps), GWASidx-1))
		w_new = np.memmap(tmpdir+'/w_tmp.dat', dtype='float128', mode='w+', shape=(len(snps)+len(new_uid), GWASidx))
		w_new[0:len(w_old),0:(GWASidx-1)] = w_old[:]
		del w_old
		w = np.memmap(tmpdir+'/w.dat', dtype='float128', mode='w+', shape=(len(snps)+len(new_uid), GWASidx))
		w[:] = w_new[:]
		del w_new
		v = np.memmap(tmpdir+'/v.dat', dtype='float128', mode='r+', shape=(len(snps)+len(new_uid), 3), order='C')
		# w = np.c_[w, [0]*len(w)]
		if len(new_uid) > 0:
			print "Detected additional "+str(len(new_uid))+" SNPs"
			tmp_idx = ArrayIn(gwas[:,0], new_uid)
			tmp = np.c_[new_uid, gwas[tmp_idx, 3:5], ["?"*(GWASidx-1)]*len(new_uid),gwas[tmp_idx, 8]]
			snps = np.r_[snps, tmp]
			n = snps[:,0].argsort()
			w[:] = w[n,:]
			v[:] = v[n,:]
			snps = snps.take(n,0)
		w[ArrayIn(snps[:,0], gwas[:,0]),GWASidx-1] = gwas[:,7]
		print "Updating variables..."
		n = ArrayIn(snps[:,0],gwas[:,0])
		v[n,0] = np.add(v[n,0], np.multiply(gwas[:,7], gwas[:,5]))
		v[n,1] = np.add(v[n,1], np.square(gwas[:,7]))
		for i in range(1,GWASidx):
			v[n,2] = np.add(v[n,2], np.multiply(w[n,i-1], w[n,GWASidx-1])*C[GWASidx-2][i-1])

		w.flush()
		v.flush()
		### align allele and add direction
		print "Aligning direction..."
		# gwas[:,6] = [x[2] if x[0]==x[1] else -1*x[2] for x in np.c_[gwas[:,3], snps[n,1], gwas[:,6]]]
		snps[n,3] = [x[0]+"+" if x[1]>0 else x[0]+"-" for x in np.c_[snps[n,3], gwas[:,6]]]
		nidx = ArrayNotIn(snps[:,0], gwas[:,0])
		snps[nidx,3] = [x+"?" for x in snps[nidx,3]]

		### check rsID
		# if "NA" is stored, take the new rsID
		# if rsID doesn't match with store one, replace with uniqID
		print "Checking rsID..."
		snps[n,4] = map(match_rsID, np.c_[gwas[:,0], gwas[:,8], snps[n,4]])

	del v,w
	return snps

### procedd GWAS
def processGWAS(config, twoside, C, chromfilt):
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

	with open(config, 'r') as inconfig:
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
				snps = processFile(gwasfile, C, snps, GWASidx, chrom, pos, a1, a2, p, effect, oddsratio, N, weight, rsID, delim, twoside, chromfilt)
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
		return GWASidx, snps

##### compute z from stored variables and combert to P #####
def computeZ(snps, nGWAS, twoside):
	v = np.memmap(tmpdir+'/v.dat', dtype='float128', mode='r+', shape=(len(snps), 3))
	w = np.memmap(tmpdir+'/w.dat', dtype='float128', mode='r+', shape=(len(snps), nGWAS))
	z = np.divide(v[:,0].astype(float), np.sqrt(np.add(v[:,1].astype(float), 2*v[:,2].astype(float))))
	if twoside:
		p = st.norm.cdf(list(-1.0*np.absolute(z)))*2
	else:
		p = st.norm.cdf(list(-1.0*z))
	Neff = np.divide(np.square(v[:,1].astype(float)), np.add(v[:,1].astype(float), 2*v[:,2].astype(float)))
	chrom = [int(re.match(r"(\d+):.+", x).group(1)) for x in snps[:,0]]
	pos = [int(re.match(r"\d+:(\d+):.+", x).group(1)) for x in snps[:,0]]

	### return chr, pos, a1, a2, rsID, z, p, weight, direction
	return np.c_[chrom, pos, snps[:,[1,2]], snps[:,4], z, p, [int(round(sum(np.square(x)))) for x in w], Neff, snps[:,3]]

def main(args):
	start_time = time.time()

	print HEADMSS

	### check arguments
	if args.config is None:
		parser.print_help()
		sys.exit("\nERROR: Config file is required.")
	if args.intercept is None:
		parser.print_help()
		sys.exit("\nERROR: Intercept file is required.")

	C = getIntercept(args.intercept)

	nGWAS, snps = processGWAS(args.config, args.twoside, C, args.chrom)

	results = computeZ(snps, nGWAS, args.twoside)
	results = results[np.lexsort((results[:,1], results[:,0]))]

	with open(args.out+".txt", 'w') as o:
		o.write("\t".join(["chr", "pos", "a1", "a2", "rsID", "z", "p", "weight", "Neff", "direction"])+"\n")
	with open(args.out+".txt", 'a') as o:
		np.savetxt(o, results, delimiter="\t", fmt="%s")
	os.system("rm -r "+tmpdir)
	print "Program run time: "+str(time.time()-start_time)

if __name__ == "__main__": main(parser.parse_args())
