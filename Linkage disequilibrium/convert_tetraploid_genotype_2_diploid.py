#!/usr/bin/env python3

import sys

i1 = sys.stdin  # .vcf.gz

for line in i1:
	if '#' in line:
		print(line.strip())
	else:
		line = line.strip().split()
		for n in range(0, 9):
			print(line[n], end='\t')
		for n in range(9, len(line)):
			m = line[n].split(":")
			if n != (len(line) - 1):
				for k in range(0, len(m)):
					if k == 0:
						if './.' in m[k]:
							print('./.', end=':')
						elif m[k].count('0') == 4:
							print('0/0', end=':')
						elif m[k].count('1') == 4:
							print('1/1', end=':')
						else:
							print('0/1', end=':')
					elif k != (len(m) - 1) and k != 0:
						print(m[k], end=':')
					else:
						print(m[k], end='\t')
			else:
				for k in range(0, len(m)):
					if k == 0:
						if './.' in m[k]:
							print('./.', end=':')
						elif m[k].count('0') == 4:
							print('0/0', end=':')
						elif m[k].count('1') == 4:
							print('1/1', end=':')
						else:
							print('0/1', end=':')
					elif k != (len(m) - 1) and k != 0:
						print(m[k], end=':')
					else:
						print(m[k])

