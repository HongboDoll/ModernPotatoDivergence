#!/usr/bin/env python3

import sys

i1 = open(sys.argv[1])  # ${vcf}.MODERATE_HIGH_gene.xls
i2 = sys.stdin  # cat 
i3 = open(sys.argv[2])  # 166_potato_starch_vs_european_Fst_region.xls
chro = str(sys.argv[3]) # chromosome, chr01

fst={}
for line in i2:
	if 'pop1' not in line:
		line = line.strip().split()
		fst[line[0]] = line[1]

region = {}
for line in i3:
	line = line.strip().split()
	if line[0] == chro:
		region[line[0]+'~'+line[1]+'~'+line[2]] = ''

for line in i1:
	if '#' not in line:
		line = line.strip().split()
		f = 0
		for k in region:
			start = int(k.split('~')[1])
			end = int(k.split('~')[2])
			if int(line[1]) >= start and int(line[1]) <= end:
				f += 1
		if f:
			diverged = 'InHighlyDivergedRegion'
		else:
			diverged = 'NotInHighlyDivergedRegion'
		if line[1] in fst:
			print('\t'.join(line), fst[line[1]], diverged, sep='\t')

