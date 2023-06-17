#!/usr/bin/env python3

import sys

i1 = open(sys.argv[1])  # 166_potato_Allchr.gatk.filter.SNP.filterMissMAF_Allchr_window_count
i2 = open(sys.argv[2])  # 166_potato_3pop_nucleotide_diversity.xls / 166_potato_starch_vs_others_nucleotide_diversity.xls
cutoff = int(sys.argv[3]) # lowest cutoff of number of SNPs; windows with lower than this value will be discarded

d = {}
for line in i1:
	line = line.strip().split()
	d[line[0]+'~'+line[1]+'~'line[2]] = int(line[3])

for line in i2:
	if 'Pos' in line:
		print(line.strip())
	else:
		line = line.strip().split()
		if line[0]+'~'+line[1]+'~'line[2] in d and d[line[0]+'~'+line[1]+'~'line[2]] < cutoff:
			print('\t'.join(line[0:3]), end='\t')
			for n in range(3, len(line)):
				if n != (len(line) - 1):
					print(0, end='\t')
				else:
					print(0)
		else:
			print('\t'.join(line))

