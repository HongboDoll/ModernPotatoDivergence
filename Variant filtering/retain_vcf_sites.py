#!/usr/bin/env python3

import sys

i1 = open(sys.argv[1])  # vcf
i2 = open(sys.argv[2])  # filter_retained_sites_chr01

d = {}
for line in i2:
    line = line.strip().split()
    d[line[0]+'~'+line[1]] = ''

for line in i1:
    if '#' in line:
        print(line.strip())
    else:
        line = line.strip().split()
        if line[0]+'~'+line[1] in d:
            print('\t'.join(line))

