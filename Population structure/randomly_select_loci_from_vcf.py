#!/usr/bin/env python3

import sys
import random

i1 = sys.stdin   # 
num = int(sys.argv[1]) # number of loci
ran = int(sys.argv[2]) # number of randomly selected loci

random_num = random.sample(range(1,num),ran)

d = {}
n = 1
for line in i1:
	if '#CHROM' in line:
		print(line.strip())
	else:
		d[n] = line.strip()
		n += 1

for k in d:
	if k in random_num:
		print(d[k])
