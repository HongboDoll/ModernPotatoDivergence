#!/usr/bin/env python3

import sys

i1 = open(sys.argv[1])  # reported_qtl_gwas_interval.bed
i3 = open(sys.argv[3])  # 166_potato_european_vs_starch_nucleotide_diversity_region_overlap_gene_information.xls
i2 = open(sys.argv[2])  # DM_v6.1.tpm.txt

qtl = {}
for line in i1:
	line = line.strip().split()
	qtl[line[3]] = [line[0], int(line[1]), int(line[2])]

tpm = {}
for line in i2:
	if 'ID' not in line:
		line = line.strip().split()
		tpm[line[0]] = line[1:]

print('Chr\tStart\tEnd\tGene\tQTL\tGO\tIPR\tSwissprot\tAra\tKEGG\tRoots_aver\tShoots\tLeaves\tStolons\tTubers\tFlowers_aver\tFruits_aver')
for line in i3:
	line = line.strip().split('\t')
	qtll = ''
	for k in qtl:
		if line[0] == qtl[k][0] and int(line[1]) >= qtl[k][1] and int(line[1]) <= qtl[k][2]:
			qtll += (k+',')
	print('\t'.join(line[0:4]), end='\t')
	if qtll:
		print(qtll, end='\t')
	else:
		print('NA', end='\t')
	print('\t'.join(line[4:]), end='\t')
	if line[3] in tpm:
		print('\t'.join(tpm[line[3]]))

