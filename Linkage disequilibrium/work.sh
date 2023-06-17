#!/bin/bash

cd /media/scratchpad_02/li322/work/04_137_atlas_population_analyses/ld/diploid_genotype_LD

#### sample1 sample2 sample... is sample list file including subgroups accessions [one column]
vcf=166_potato_Allchr.gatk.filter.SNP.filterMissMAF.vcf.gz
diploid_vcf=166_potato_Allchr.gatk.filter.SNP.filterMissMAF.diploid_genotype.vcf

zcat 166_potato_Allchr.gatk.filter.SNP.filterMissMAF.vcf.gz | ./convert_tetraploid_genotype_2_diploid.py > ${diploid_vcf}

PopLDdecay -InVCF ${diploid_vcf} -SubPop pop1_starch.xls -MaxDist 500 -MAF 0.05 -Miss 0.1 -OutStat ${diploid_vcf}.pop1_starch.stat.gz
PopLDdecay -InVCF ${diploid_vcf} -SubPop pop2_european.xls -MaxDist 500 -MAF 0.05 -Miss 0.1 -OutStat ${diploid_vcf}.pop2_european.stat.gz
PopLDdecay -InVCF ${diploid_vcf} -SubPop pop3_american.xls -MaxDist 500 -MAF 0.05 -Miss 0.1 -OutStat ${diploid_vcf}.pop3_american.stat.gz

PopLDdecay -InVCF ${diploid_vcf} -MaxDist 500 -MAF 0.05 -Miss 0.1 -OutStat ${diploid_vcf}.all.stat.gz

echo -e "${diploid_vcf}.pop1_starch.stat.gz\tpop1_starch\n${diploid_vcf}.pop2_european.stat.gz\tpop2_european\n${diploid_vcf}.pop3_american.stat.gz\tpop3_american" > 166_potato_filter.LDdecayResult.list
perl /media/bulk_01/users/li322/software/PopLDdecay-3.42/bin/Plot_MultiPop.pl -inList  166_potato_filter.LDdecayResult.list -keepR -output 166_potato_filter_population_LDdecay
perl /media/bulk_01/users/li322/software/PopLDdecay-3.42/bin/Plot_OnePop.pl -inFile ${diploid_vcf}.all.stat.gz -keepR -output 166_potato_filter_all_LDdecay

./plot_LD_popLDdecay.R 166_potato_filter_population_LDdecay.pop1_starch 166_potato_filter_population_LDdecay.pop2_european 166_potato_filter_population_LDdecay.pop3_american 166_potato_filter_population_LDdecay.pdf
