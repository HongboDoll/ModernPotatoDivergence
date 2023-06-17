#!/usr/bin/env python3

import sys

i1 = open(sys.argv[1])  # 

for line in i1:
	if '#' not in line:
		line = line.strip().split()
		print('\t'.join(line[0:7]), 'DP=2000', 'GT', sep='\t', end='\t')
		for n in range(9, len(line)):
#			print(line[n].split(':')[0], end='\t')
			if '.' in line[n].split(':')[0]:
				print('./.', end='\t')
			elif '.' not in line[n].split(':')[0] and line[n].split(':')[0] != '0/0/0/0' and line[n].split(':')[0] != '1/1/1/1/':
				print('0/1', end='\t')
			else:
				print(line[n].split(':')[0][0:3], end='\t')
		print('')
	else:
		print(line.strip())

