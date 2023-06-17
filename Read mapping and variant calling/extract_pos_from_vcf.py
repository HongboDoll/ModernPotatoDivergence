#!/usr/bin/env python3

import sys

i1 = sys.stdin  # zcat 166_potato_Allchr.gatk.filter.indel.vcf.snpeff.gz
i2 = open(sys.argv[1])  # 166_potato_Allchr.gatk.filter.${i}.${n}.private.xls

pos = {}
for line in i2:
	line = line.strip().split()
	pos[line[0]+'~'+line[1]] = ''

for line in i1:
	if '#' not in line:
		line  = line.strip().split()
		if line[0]+'~'+line[1] in pos:
			print('\t'.join(line))

