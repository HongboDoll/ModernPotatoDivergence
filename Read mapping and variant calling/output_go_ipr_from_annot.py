#!/usr/bin/env python3

import sys

i1 = open(sys.argv[1])  # 
o1 = open(sys.argv[2], 'w')  # go
o2 = open(sys.argv[3], 'w')  # ipr

for line in i1:
	if '#Chr' not in line:
		line = line.strip().split('\t')
		go = line[4]
		ipr = line[5]
		print(line[3], end='\t')
		o1.write('%s\t' % line[3])
		for k in range(0, len(go.split(';')[:-1])):
			if k != (len(go.split(';')[:-1]) - 1):
				print(go.split(';')[:-1][k].split()[0], end='\t')
				o1.write('%s,' % go.split(';')[:-1][k].split()[0])
			else:
				print(go.split(';')[:-1][k].split()[0])
				o1.write('%s\n' % go.split(';')[:-1][k].split()[0])
		print()
		o1.write('\n')
		o2.write('%s\t' % line[3])
		for k in range(0, len(ipr.split(';')[:-1])):
			if k != (len(ipr.split(';')[:-1]) - 1):
				o2.write('%s\t' % ipr.split(';')[:-1][k].split()[0])
			else:
				o2.write('%s\n' % ipr.split(';')[:-1][k].split()[0])
		o2.write('\n')
