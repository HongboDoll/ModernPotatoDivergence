#!/usr/bin/env python3

import sys

i1 = open(sys.argv[1])  # vcf
i2 = float(sys.argv[2]) # proportion of samples being heterozygous, 0.5, 0.8

for line in i1:
    if '#' not in line:
        line = line.strip().split()
        alt = 0
        ref = 0
        mis = 0
        het = 0
        for n in range(9, len(line)):
            gt = line[n].split(':')[0]
            if len(gt.split('/')) == 2:
                if gt == '0/0':
                    ref += 1
                elif gt == '1/1':
                    alt += 1
                elif gt == './.':
                    mis += 1
                else:
                    het += 1
            elif len(gt.split('/')) == 4:
                if gt == '0/0/0/0':
                    ref += 1
                elif gt == '1/1/1/1':
                    alt += 1
                elif gt == './././.':
                    mis += 1
                else:
                    het += 1
                    
        if alt+ref+het != 0:
            if float(het/float(alt+het)) > i2:
                print('\t'.join(line[0:3]))

