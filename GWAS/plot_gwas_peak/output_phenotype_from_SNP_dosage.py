#!/usr/bin/env python3

import sys

i1 = open(sys.argv[1])  # leading_snp.vcf
i2 = open(sys.argv[2])  # 166_potato_50_phenotype_population_structure.xls

snp_dosage = {}
species = {}
for line in i1:
	if '#CHROM' in line:
		line = line.strip().split()
		for k in range(9, len(line)):
			species[k] = line[k]
	else:
		if '#' not in line:
			line = line.strip().split()
			for k in range(9, len(line)):
				gt = line[k].split(':')[0]
				if '.' not in gt:
					if gt not in snp_dosage:
						snp_dosage[gt] = []
						snp_dosage[gt].append(species[k])
					else:
						snp_dosage[gt].append(species[k])

pheno = {}
for line in i2:
	if 'Name' not in line:
		line = line.strip().split()
		pheno[line[0]] = line[1]


snp_dosage_pheno = {}
for k in snp_dosage:
	snp_dosage_pheno[k] = []
	for j in snp_dosage[k]:
		if j in pheno and pheno[j] != "NA":
			snp_dosage_pheno[k].append(pheno[j])

for k in sorted(snp_dosage_pheno):
	print(k, end='\t')
	for j in snp_dosage_pheno[k]:
		print(j, end='\t')
	print('')

