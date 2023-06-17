#!/usr/bin/env python3

import sys

i1 = sys.stdin  # genotype/115_v3_chr1-7_snp_ID_filter_LDpruned_10.genotype

d = {'A':'1', 'T':'2', 'C':'3', 'G':'4'}
for line in i1:
	line = line.strip().split()
	if 'REF' in line[0]:
		print(line[2], end=' ')
	else:
		if '/' in line[2]:
			if len(line[2].split('/')) == 2:
				for k in line[2].split(':')[0].split('/'):
					if k == '0':
						print(d[line[0]]+' '+d[line[0]], end=' ')
					elif k == '1':
						print(d[line[1]]+' '+d[line[1]], end=' ')
					elif k == '.':
						print('-9 -9', end=' ')
						
			elif len(line[2].split('/')) == 4:
				for k in line[2].split(':')[0].split('/'):
					if k == '0':
						print(d[line[0]], end=' ')
					elif k == '1':
						print(d[line[1]], end=' ')
					elif k == '.':
						print('-9', end=' ')

		elif '|' in line[2]:
			if len(line[2].split('|')) == 2:
				for k in line[2].split(':')[0].split('|'):
					if k == '0':
						print(d[line[0]]+' '+d[line[0]], end=' ')
					elif k == '1':
						print(d[line[1]]+' '+d[line[1]], end=' ')
					elif k == '.':
						print('-9 -9', end=' ')

			elif len(line[2].split('|')) == 4:
				for k in line[2].split(':')[0].split('|'):
					if k == '0':
						print(d[line[0]], end=' ')
					elif k == '1':
						print(d[line[1]], end=' ')
					elif k == '.':
						print('-9', end=' ')
print('')
