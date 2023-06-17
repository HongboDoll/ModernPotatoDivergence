#!/usr/bin/env python3

import sys

i1 = sys.stdin  # 

for line in i1:
	if '#' in line:
		print(line.strip())
	else:
		line = line.strip().split()
		print('\t'.join(line[:8]), 'GT', sep='\t', end='\t')
		for k in range(9, len(line)):
			if k != (len(line) - 1):
				print(line[k].split(':')[0], end='\t')
			else:
				print(line[k].split(':')[0])

