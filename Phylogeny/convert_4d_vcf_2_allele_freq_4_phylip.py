#!/usr/bin/env python3

import sys

i1 = open(sys.argv[1])  # 115_v3_chr1-7_snp_ID_filter.vcf.4d
i2 = open(sys.argv[2]) # species name correspondance to meet the ten characters requirement

n_sample = 0
n_loci = 0
name = []
d = {}
for line in i1:
	if '#CHROM' in line:
		n_sample = len(line.strip().split()[9:])
		for k in line.strip().split()[9:]:
			name.append(k)
	elif '#CHROM' not in line and '#' not in line:
		n_loci += 1
		line = line.strip().split()
		for k in range(0, len(line[9:])):
			if name[k] not in d:
				d[name[k]] = []
				d[name[k]].append('%.4f' % (float(line[9:][k].split(':')[0].split('/').count('0'))/len(line[9:][k].split(':')[0].split('/'))))
			else:
				d[name[k]].append('%.4f' % (float(line[9:][k].split(':')[0].split('/').count('0'))/len(line[9:][k].split(':')[0].split('/'))))
			
corres = {}
for line in i2:
	line = line.strip().split()
	corres[line[1]] = line[0]

print('    '+str(n_sample)+'    '+str(n_loci))
for k in range(0, n_loci):
	print('2', end=' ')
print('')

for k in d:
	print(corres[k], ' '.join(d[k]), sep='          ')
	

