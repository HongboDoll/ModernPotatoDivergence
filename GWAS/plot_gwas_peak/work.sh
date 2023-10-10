#!/bin/bash

cd /media/scratchpad_02/li322/work/07_gwas/plot_blue_uww_maturity

#./manhanttan_CMplot.R plot/166_potato_50_phenotype_GWASpoly_results.xls_Allchr_Blue_uww_additive_4_plot Blue_uww_additive_gwas_plot.pdf 1.50E-7

#./manhanttan_CMplot.R plot/166_potato_50_phenotype_GWASpoly_results.xls_Allchr_Maturity_additive_4_plot Maturity_additive_gwas_plot.pdf 1.50E-7

#rm blue_gene_coordinate.xls
#cat blue_gene.xls | while read i
#do
#	id=`grep $i DM_v6.1.gff3 | awk '$3=="mRNA"' | awk '{print $NF}' | awk -F ';' '{print $1}' | awk -F '=' '{print $2}'`
#	grep $i DM_v6.1.gff3 | awk '$3=="mRNA"' | awk 'BEGIN{OFS="\t"}{print "'"$id"'",$1,$4,$5,$7,"1","blue"}' >> blue_gene_coordinate.xls
#done
#
#rm maturity_gene_chr05_coordinate.xls
#cat maturity_gene_chr05.xls | while read i
#do
#    id=`grep $i DM_v6.1.gff3 | awk '$3=="mRNA"' | awk '{print $NF}' | awk -F ';' '{print $1}' | awk -F '=' '{print $2}'`
#    grep $i DM_v6.1.gff3 | awk '$3=="mRNA"' | awk 'BEGIN{OFS="\t"}{print "'"$id"'",$1,$4,$5,$7,"1","blue"}' >> maturity_gene_chr05_coordinate.xls
#done
#
#
#rm maturity_gene_chr09_coordinate.xls
#cat maturity_gene_chr09.xls | while read i
#do
#    id=`grep $i DM_v6.1.gff3 | awk '$3=="mRNA"' | awk '{print $NF}' | awk -F ';' '{print $1}' | awk -F '=' '{print $2}'`
#    grep $i DM_v6.1.gff3 | awk '$3=="mRNA"' | awk 'BEGIN{OFS="\t"}{print "'"$id"'",$1,$4,$5,$7,"1","blue"}' >> maturity_gene_chr09_coordinate.xls
#done
#
#
#cat <(grep '#' leading_snp.vcf) <(awk '$3=="SNP.chr06_50213402"' leading_snp.vcf) > tmp.vcf
#cut -f '1,46' 166_potato_50_phenotype_population_structure.xls > tmp_pheno.xls
#
#./output_phenotype_from_SNP_dosage.py tmp.vcf tmp_pheno.xls > blue_leading_snp_dosage_phenotype.xls
#t.R blue_leading_snp_dosage_phenotype.xls blue_leading_snp_dosage_phenotype.xls_t.xls
#
#cat <(grep '#' leading_snp.vcf) <(awk '$3=="SNP.chr05_4810976"' leading_snp.vcf) > tmp.vcf
#cut -f '1,4' 166_potato_50_phenotype_population_structure.xls > tmp_pheno.xls
#
#./output_phenotype_from_SNP_dosage.py tmp.vcf tmp_pheno.xls > maturity_chr05_leading_snp_dosage_phenotype.xls
#t.R maturity_chr05_leading_snp_dosage_phenotype.xls maturity_chr05_leading_snp_dosage_phenotype.xls_t.xls
#
#cat <(grep '#' leading_snp.vcf) <(awk '$3=="SNP.chr09_59264287"' leading_snp.vcf) > tmp.vcf
#cut -f '1,4' 166_potato_50_phenotype_population_structure.xls > tmp_pheno.xls
#
#./output_phenotype_from_SNP_dosage.py tmp.vcf tmp_pheno.xls > maturity_chr09_leading_snp_dosage_phenotype.xls
#t.R maturity_chr09_leading_snp_dosage_phenotype.xls maturity_chr09_leading_snp_dosage_phenotype.xls_t.xls
#
#rm blue_leading_snp_dosage_phenotype.xls_4_anova
#cat blue_leading_snp_dosage_phenotype.xls | while read i
#do
#	n=`echo $i | awk '{print $1}'`
#	echo $i | sed 's/ /\n/g' | awk 'NR!=1' | while read s
#	do
#		echo -e "${n}\t${s}" >> blue_leading_snp_dosage_phenotype.xls_4_anova
#	done
#done
#
#
#rm maturity_chr05_leading_snp_dosage_phenotype.xls_4_anova
#cat maturity_chr05_leading_snp_dosage_phenotype.xls | while read i
#do
#    n=`echo $i | awk '{print $1}'`
#    echo $i | sed 's/ /\n/g' | awk 'NR!=1' | while read s
#    do
#        echo -e "${n}\t${s}" >> maturity_chr05_leading_snp_dosage_phenotype.xls_4_anova
#    done
#done
#
#rm maturity_chr09_leading_snp_dosage_phenotype.xls_4_anova
#cat maturity_chr09_leading_snp_dosage_phenotype.xls | while read i
#do
#    n=`echo $i | awk '{print $1}'`
#    echo $i | sed 's/ /\n/g' | awk 'NR!=1' | while read s
#    do
#        echo -e "${n}\t${s}" >> maturity_chr09_leading_snp_dosage_phenotype.xls_4_anova
#    done
#done

# phylogeny of EID1 homologs among tomato potato and ara
cat potato_EID1_homologs.fa Ara_EID1_homologs.fa tomato_EID1_homologs.fa > potato_tomato_ARA_EID1_homologs.fa

mafft --thread 10 --auto potato_tomato_ARA_EID1_homologs.fa > potato_tomato_ARA_EID1_homologs.fa.out

fa2phy.py -i potato_tomato_ARA_EID1_homologs.fa.out -o potato_tomato_ARA_EID1_homologs.fa.phylip

rm potato_tomato_ARA_EID1_homologs.fa.phylip.* RAxML_*
#iqtree -s potato_tomato_ARA_EID1_homologs.fa.phylip -m MFP -b 100 -safe -nt 10
raxmlHPC-PTHREADS-SSE3 -f a -x 5 -p 5 -# 100 -m PROTGAMMAJTT -s potato_tomato_ARA_EID1_homologs.fa.phylip -n potato_tomato_ARA_EID1_homologs.fa -T 10

