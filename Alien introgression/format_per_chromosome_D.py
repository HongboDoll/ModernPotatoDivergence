#!/usr/bin/env python3

import sys

i1 = open(sys.argv[1])  # Dtrios/output_chr${i}_BBAA.txt
i2 = open(sys.argv[2])  # whole_genome_combined_D_Z_order.xls
chrom = sys.argv[3] # chr1

spe = []
for line in i2:
	line = line.strip().split()
	spe.append(line[1])

d = {}
for line in i1:
	line = line.strip().split()
	if line[0] == 'European' and line[1] == 'Starch':
		d[line[2]] = [line[3], line[4], line[6]]
print(chrom, end='\t')
for n in range(0, len(spe)):
	if n != (len(spe) - 1):
		if spe[n] in d:
			print(d[spe[n]][0], end='\t')
		else:
			print('NA', end='\t')
	else:
		if spe[n] in d:
			print(d[spe[n]][0])
		else:
			print('NA')

