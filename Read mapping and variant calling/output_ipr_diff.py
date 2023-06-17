#!/usr/bin/env python3

import sys

i1 = open(sys.argv[1])  #  core_gene_IPR.xls
i2 = open(sys.argv[2])  #  dispensable_gene_IPR.xls
i3 = open(sys.argv[3])  #  Csa37_80_64_9930_IPR_infor.xls

ipr = {}
for line in i3:
	line = line.strip()
	sline = line.split()
	ipr[sline[0]] = line[10:]
	
core_ipr = {}
for line in i1:
	line = line.strip().split()
	core_ipr[line[0]] = line[1:]

dis_ipr = {}
for line in i2:
	line = line.strip().split()
	dis_ipr[line[0]] = line[1:]
ipr_core_num = {}
ipr_dis_num = {}
for k in ipr:
	n = 0
	for i in core_ipr:
		if k in core_ipr[i]:
			n += 1
	ipr_core_num[k] = n
	n = 0
	for i in dis_ipr:
		if k in dis_ipr[i]:
			n += 1
	ipr_dis_num[k] = n

print('#core_num\tcore_ratio\tdis_num\tdis_ratio\tdescription\tnum_core_ipr\tnum_dis_ipr')
for k in ipr_core_num:
	if k in ipr_dis_num:
		print(ipr_core_num[k],'%.4f' % (ipr_core_num[k]/len(core_ipr)*100), ipr_dis_num[k],'%.4f' % (ipr_dis_num[k]/len(dis_ipr)*100), k+':'+ipr[k],len(core_ipr),len(dis_ipr),sep='\t')
		

