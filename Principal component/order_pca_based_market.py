#!/usr/bin/env python3

import sys

i1 = open(sys.argv[1])  # 166_potato_market_order.xls
i2 = open(sys.argv[2])  # 166_potato_Allchr.gatk.filter.SNP.filterMissMAF.4d.pca

market = {}
for line in i1:
    line = line.strip().split('~')
    if line[1] not in market:
        market[line[1]] = [line[0]]
    else:
        market[line[1]].append(line[0])

pca = {}
for line in i2:
    line = line.strip().split()
    pca[line[0]] = line

for k in market:
    for j in market[k]:
        if j in pca:
            print('\t'.join(pca[j]))

