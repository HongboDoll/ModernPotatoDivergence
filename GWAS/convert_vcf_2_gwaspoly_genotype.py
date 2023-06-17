#!/usr/bin/env python3

import sys

i1 = open(sys.argv[1])  # 166_potato_Allchr.gatk.filter.SNP.filterMissMAF.vcf
#i2 = int(sys.argv[2])  # 4, ploidy level

for line in i1:
	if '#CHROM' in line:
		line = line.strip().split()
		print('Marker\tChrom\tPosition', '\t'.join(line[9:]), sep='\t')
	elif '#CHROM' not in line and '#' not in line:
		line = line.strip().split()
		print(line[2], line[0], line[1], sep='\t', end='\t')
		geno = {'0':line[3], '1':line[4]}
		for n in range(9, len(line)):
			gt = line[n].split(':')[0]
			if n != (len(line) - 1):
				if '.' not in gt:
					for k in gt.split('/'):
						print(geno[k], end='')
					print('', end='\t')
				else:
					print('NA', end='\t')
			else:
				if '.' not in gt:
					for k in gt.split('/'):
						print(geno[k], end='')
					print('')
				else:
					print('NA')

