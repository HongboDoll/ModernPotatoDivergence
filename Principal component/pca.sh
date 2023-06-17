#!/bin/bash

cd /media/scratchpad_02/li322/work/04_137_atlas_population_analyses/pca

#cd ..
#cat head_vcf <(grep -v '#' 166_potato_Allchr.gatk.filter.SNP.filterMissMAF.vcf | sort -k1,1 -k2,2n) > t && mv t 166_potato_Allchr.gatk.filter.SNP.filterMissMAF.vcf
#cat head_vcf <(grep -v '#' 166_potato_Allchr.gatk.filter.SNP.filterMissMAF.4d.vcf | sort -k1,1 -k2,2n) > t && mv t 166_potato_Allchr.gatk.filter.SNP.filterMissMAF.4d.vcf
#cat <(grep '#' 166_potato_Allchr.gatk.filter.SNP.filterMissMAF.4d.vcf) <(paste <(grep -v '#' 166_potato_Allchr.gatk.filter.SNP.filterMissMAF.4d.vcf | cut -f '1-2') <(grep -v '#' 166_potato_Allchr.gatk.filter.SNP.filterMissMAF.4d.vcf | awk 'BEGIN{n=1}{print "SNP"n;n+=1}') <(grep -v '#' 166_potato_Allchr.gatk.filter.SNP.filterMissMAF.4d.vcf | cut -f '4-999')) > t && mv t 166_potato_Allchr.gatk.filter.SNP.filterMissMAF.4d.vcf

#cat <(grep '#' 166_potato_Allchr.gatk.filter.SNP.filterMissMAF.vcf) <(paste <(grep -v '#' 166_potato_Allchr.gatk.filter.SNP.filterMissMAF.vcf | cut -f '1-2') <(grep -v '#' 166_potato_Allchr.gatk.filter.SNP.filterMissMAF.vcf | awk 'BEGIN{n=1}{print "SNP"n;n+=1}') <(grep -v '#' 166_potato_Allchr.gatk.filter.SNP.filterMissMAF.vcf | cut -f '4-999')) > t && mv t 166_potato_Allchr.gatk.filter.SNP.filterMissMAF.vcf

#cd -
#
#### 4d SNP dataset 
vcf=166_potato_Allchr.gatk.filter.SNP.filterMissMAF.4d.vcf
end=`grep -v '#' ${vcf} |head -1| awk '{print NF}'`

grep '#CHROM' head_vcf | cut -f '10-999'|sed 's/\t/\n/g'|awk '{print $1"\t4"}' > accession_ploidy.xls

cat <(grep '#CHROM' ${vcf}) <(grep -v '#' ${vcf}) > $(basename ${vcf} .vcf).genotype

rm -rf struc_tmp; mkdir struc_tmp
for n in `seq 10 $end`
    do
    cut -f 4,5,$n $(basename ${vcf} .vcf).genotype | ./convert_vcf_2_structure.py > struc_tmp/$(basename ${vcf} .vcf)_${n}.structure
done
	
cat struc_tmp/*structure | sort -k1,1 > $(basename ${vcf} .vcf).structure

cat <(echo -e " \c\n") <(grep -v '#' $vcf  | awk '{print $3}' | sed ':a;N;s/\n/ /g;ta') > head
cat head <(./convert_structrue_2_adegenet.py $(basename ${vcf} .vcf).structure accession_ploidy.xls) > $(basename ${vcf} .vcf).adegenet

./pca.R $(basename ${vcf} .vcf).adegenet $(basename ${vcf} .vcf).pca

./order_pca_based_market.py 166_potato_market_order.xls  166_potato_Allchr.gatk.filter.SNP.filterMissMAF.4d.pca > 166_potato_Allchr.gatk.filter.SNP.filterMissMAF.4d.pca_order

./PCA_plot_variance_explained.R 166_potato_Allchr.gatk.filter.SNP.filterMissMAF.4d.pca_order 166_potato_Allchr.gatk.filter.SNP.filterMissMAF.4d.pca_order.pdf
