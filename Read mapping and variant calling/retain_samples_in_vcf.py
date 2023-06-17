#!/usr/bin/env python3

import sys

i1 = open(sys.argv[1])  # merged.vcf without genotype
i2 = open(sys.argv[2])  # 12_sol_wild.xls, sampels to be retained in the resulted vcf
i3 = int(sys.argv[3]) # ploidy level 2/4

acc1 = {}
for line in i2:
	acc1[line.strip()] = ''

corres = {}
retain = {}
for line in i1:
	if '#' in line and '#CHROM' not in line:
		print(line.strip())
	elif '#' in line and '#CHROM' in line:
		line = line.strip().split()
		print('\t'.join(line[0:9]),end='\t')
		for n in range(9, len(line)):
			corres[n] = line[n]
			if corres[n] in acc1:
				retain[n] = ''
				print(line[n], end='\t')
		print('')
	else:
		line = line.strip().split()
		f = 0
		for k in retain:
			if i3 == 2:
				if line[k].split(':')[0] != '0/0' and line[k].split(':')[0] != './.':
					f += 1
			elif i3 == 4:
				if line[k].split(':')[0] != '0/0/0/0' and line[k].split(':')[0] != './././.':
					f += 1
		if f:
			print('\t'.join(line[0:9]),end='\t')
			for k in retain:
				print(line[k], end='\t')
			print('')


