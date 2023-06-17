#!/usr/bin/env python3

import sys

i1 = open(sys.argv[1])  # whole_genome_combined_BBAA.txt
i2 = open(sys.argv[2])  # wild_order.xls

order = []
for line in i2:
	order.append(line.strip())

d = {}
for line in i1:
	line = line.strip().split()
	if line[0] == 'European' and line[1] == 'Starch':
		d[line[2]] = [line[3], line[4]]

for k in order:
	if k in d:
		print('S. '+k, d[k][0], d[k][1], sep='\t')

