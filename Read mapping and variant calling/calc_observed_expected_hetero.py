#!/usr/bin/env python3

import sys

i1 = sys.stdin  # vcf

for line in i1:
	if '#' not in line:
		line = line.strip().split()
		alt = 0
		ref = 0
		mis = 0
		het = 0
		all_alleles = 0
		ref_alleles = 0
		alt_alleles = 0
		for n in range(9, len(line)):
			gt = line[n].split(':')[0].replace('|', '/')
			if len(gt.split('/')) == 2:
				if gt == '0/0':
					ref += 1
				elif gt == '1/1':
					alt += 1
				elif gt == './.':
					mis += 1
				else:
					het += 1
			elif len(gt.split('/')) == 4:
				if gt == '0/0/0/0':
					ref += 1
				elif gt == '1/1/1/1':
					alt += 1
				elif gt == './././.':
					mis += 1
				else:
					het += 1
			if gt != './.' or gt != './././.':
				for k in gt.split('/'):
					all_alleles += 1
					if k == "0":
						ref_alleles += 1
					elif k == "1":
						alt_alleles += 1
		if all_alleles != 0:
			exp_heter = 2 * (float(ref_alleles)/all_alleles) * (float(alt_alleles)/all_alleles)
		else:
			exp_heter = 0
		if alt+het+ref != 0:
			ob_heter = float(het/float(alt+het+ref))
		else:
			ob_heter = 0
		print('\t'.join(line[0:5]), ob_heter, exp_heter, sep='\t')

