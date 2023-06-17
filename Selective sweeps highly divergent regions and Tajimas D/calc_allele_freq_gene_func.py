#!/usr/bin/env python3

import sys

i1 = open(sys.argv[1])  # pop1_starch.xls
i2 = open(sys.argv[2])  # pop2_european.xls
i3 = open(sys.argv[3])  # 166_potato_Allchr.gatk.filter.SNP.filterMissMAF.vcf.MODERATE_HIGH.snpeff
i4 = open(sys.argv[4])  # 166_potato_Allchr.gatk.filter.SNP.filterMissMAF.vcf.MODERATE_HIGH_gene_fst0.3_inDivergedRegions.xls
i5 = open(sys.argv[5])  # DM_v6.1_GO_IPR_Swiss_Ara_KEGG.txt
i6 = open(sys.argv[6])  # reported_qtl_gwas_interval.bed
i7 = open(sys.argv[7])  # DM_v6.1.tpm.txt

starch = {}
for line in i1:
	starch[line.strip()] = ''

european = {}
for line in i2:
	european[line.strip()] = ''

starch_n = {}
european_n = {}
allele_freq = {}
print('Chr\tPOS\tID\tREF\tALT\tGene\tFST\tALT allele frequency in starch\tALT allele frequency in others\tGO\tIPR\tSwissprot\tAra\tKEGG\tQTL\tRoots_aver\tShoots\tLeaves\tStolons\tTubers\tFlowers_aver\tFruits_aver')
for line in i3:
	if '#CHROM' in line:
		line = line.strip().split()
		for n in range(9, len(line)):
			if line[n] in starch:
				starch_n[n] = ''
			elif line[n] in european:
				european_n[n] = ''
	elif '#CHROM' not in line and '#' not in line:
		line = line.strip().split()
		count_starch = 0
		count_others = 0
		count_all_starch = 0
		count_all_others = 0
		for n in range(9, len(line)):
			gt = line[n].split(':')[0]
			if n in starch_n:
		#		count_all_starch += 4
				if '.' not in gt:
					count_all_starch += 4
					count_starch += gt.split('/').count('1')
			elif n in european_n:
#				count_all_others += 4
				if '.' not in gt:
					count_all_others += 4
					count_others += gt.split('/').count('1')
		if count_all_starch != 0 and count_all_others != 0:
			allele_freq[line[0]+'~'+line[1]] = ['%.4f' % (count_starch/float(count_all_starch)), '%.4f' % (count_others/float(count_all_others)), '%.4f' % ((count_starch+count_others)/float(count_all_starch+count_all_others))]
		elif count_all_starch == 0 and count_all_others != 0:
			allele_freq[line[0]+'~'+line[1]] = [0, '%.4f' % (count_others/float(count_all_others)), '%.4f' % ((count_starch+count_others)/float(count_all_starch+count_all_others))]
		elif count_all_starch != 0 and count_all_others == 0:
			allele_freq[line[0]+'~'+line[1]] = ['%.4f' % (count_starch/float(count_all_starch)), 0, '%.4f' % ((count_starch+count_others)/float(count_all_starch+count_all_others))]
		elif count_all_starch == 0 and count_all_others == 0:
			allele_freq[line[0]+'~'+line[1]] = [0, 0, 0]


func = {}
for line in i5:
	if '#Chr' not in line:
		line = line.strip().split('\t')
		func[line[3]] = line[4:]

qtl = {}
for line in i6:
	line = line.strip().split()
	qtl[line[3]] = [line[0], int(line[1]), int(line[2])]

tpm = {}
for line in i7:
	if 'ID' not in line:
		line = line.strip().split()
		tpm[line[0]] = line[1:]

for line in i4:
	line = line.strip().split()
	print('\t'.join(line[0:7]), end='\t')
	if line[0]+'~'+line[1] in allele_freq:
		#p1 = float(allele_freq[line[0]+'~'+line[1]][0])
		#p2 = float(allele_freq[line[0]+'~'+line[1]][1])
		#p3 = float(allele_freq[line[0]+'~'+line[1]][2])
		#fst = (2*p3*(1-p3) - ((2*p1*(1-p1)+2*p2*(1-p2))/2))/(2*p3*(1-p3))
		print(allele_freq[line[0]+'~'+line[1]][0], allele_freq[line[0]+'~'+line[1]][1], sep='\t', end='\t')
	if line[5] in func:
		print('\t'.join(func[line[5]]), end='\t')
	qtll = ''
	for k in qtl:
		if line[0] == qtl[k][0] and int(line[1]) >= qtl[k][1] and int(line[1]) <= qtl[k][2]:
			qtll += (k+',')
	if qtll:
		print(qtll, end='\t')
	else:
		print('NA', end='\t')
	if line[5] in tpm:
		print('\t'.join(tpm[line[5]]))








