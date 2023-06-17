#!/usr/bin/env python3

import sys,re

i1 = open(sys.argv[1])  # func annotation
spe = sys.argv[2] # no for no header

for line in i1:
	if 'Gene_ID' not in line:
		sline = line.strip().split('\t')
		go = re.findall('IPR\d{6}',sline[5])
		if spe != 'no':
			print(spe+'_'+sline[3],end='\t')
			for i in go:
				print(i,end='\t')
			print('')
		else:
			print(sline[3],end='\t')
			for i in go:
				print(i,end='\t')
			print('')

