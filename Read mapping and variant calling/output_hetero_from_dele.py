#!/usr/bin/env python3

import sys

i1 = open(sys.argv[1])  # 
i2 = sys.stdin  # 

d = {}
for line in i1:
    line = line.strip().split()
    d[line[0]+'~'+line[1]+'~'+line[2]+'~'+line[3]] = [line[4], line[8]]
    
for line in i2:
    line = line.strip().split()
    if line[0]+'~'+line[1]+'~'+line[3]+'~'+line[4] in d:
        print(line[0]+'\t'+line[1]+'\t'+line[3]+'\t'+line[4]+'\t'+line[5]+'\t'+line[6], d[line[0]+'~'+line[1]+'~'+line[3]+'~'+line[4]][0], d[line[0]+'~'+line[1]+'~'+line[3]+'~'+line[4]][1], sep='\t')

