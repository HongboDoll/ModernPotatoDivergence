#!/usr/bin/env python3

import sys

i1 = open(sys.argv[1])  # 115_v3_chr1-7_snp_ID_filter_LDpruned.structure
i2 = open(sys.argv[2])  # accession_ploidy.xls

d = {}
for line in i2:
	line = line.strip().split()
	d[line[0]] = int(line[1])
	
n = 1
for line in i1:
	line = line.strip().split()
	if n == 1:
		print(line[0], end=' ')
		if d[line[0]] == 2:
			for k in range(1, len(line[1:]), 4):
				if line[k] == '-9':
					print('NA', end=' ')
				else:
					print(line[k]+'/'+line[k+2], end=' ')
		elif d[line[0]] == 4:
			for k in range(1, len(line[1:]), 4):
				if line[k] == '-9':
					print('NA', end=' ')
				else:
					print(line[k]+'/'+line[k+1]+'/'+line[k+2]+'/'+line[k+3], end=' ')

		n += 1
		print('')
	else:
		print(line[0], end=' ')
		if d[line[0]] == 2:
			for k in range(1, len(line[1:]), 4):
				if line[k] == '-9':
					print('NA', end=' ')
				else:
					print(line[k]+'/'+line[k+2], end=' ')
		elif d[line[0]] == 4:
			for k in range(1, len(line[1:]), 4):
				if line[k] == '-9':
					print('NA', end=' ')
				else:
					print(line[k]+'/'+line[k+1]+'/'+line[k+2]+'/'+line[k+3], end=' ')

		n += 1
		print('')

