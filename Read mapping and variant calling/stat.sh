#!/bin/bash

cd /media/scratchpad_02/li322/work/02_137_atlas_potato_reseq/stat

echo "all SNPs"
grep -v '#' 166_potato_Allchr.gatk.filter.SNP.vcf | awk '{print $1"\t"$2}' | sort -k1,1 -k2,2n > 166_potato_Allchr.gatk.filter.SNP.xls
wc -l 166_potato_Allchr.gatk.filter.SNP.xls

echo "all InDels"
grep -v '#' 166_potato_Allchr.gatk.filter.indel.vcf | awk '{print $1"\t"$2}' | sort -k1,1 -k2,2n > 166_potato_Allchr.gatk.filter.indel.xls
wc -l 166_potato_Allchr.gatk.filter.indel.xls


for i in pop1_starch pop2_european pop3_american
do
echo -e  "${i}\tSNP"
./retain_samples_in_vcf.py 166_potato_Allchr.gatk.filter.SNP.vcf ${i}.xls 4 |grep -v '#' | awk '{print $1"\t"$2}' | sort -k1,1 -k2,2n > 166_potato_Allchr.gatk.filter.SNP.${i}.xls
wc -l 166_potato_Allchr.gatk.filter.SNP.${i}.xls

echo -e "${i}\tInDel"
./retain_samples_in_vcf.py 166_potato_Allchr.gatk.filter.indel.vcf ${i}.xls 4 |grep -v '#' | awk '{print $1"\t"$2}' | sort -k1,1 -k2,2n > 166_potato_Allchr.gatk.filter.indel.${i}.xls
wc -l 166_potato_Allchr.gatk.filter.indel.${i}.xls

done

#### group specific variants
cat pop2_european.xls pop3_american.xls > pop23.xls
cat pop1_starch.xls pop2_european.xls > pop12.xls
cat pop1_starch.xls pop3_american.xls > pop13.xls

for i in SNP indel
do
echo -e "starch\t${i} private"
./retain_samples_in_vcf.py 166_potato_Allchr.gatk.filter.${i}.vcf pop23.xls 4|grep -v '#' | awk '{print $1"\t"$2}' | sort -k1,1 -k2,2n > 166_potato_Allchr.gatk.filter.${i}.pop23.xls

echo -e "european\t${i} private"
./retain_samples_in_vcf.py 166_potato_Allchr.gatk.filter.${i}.vcf pop13.xls 4|grep -v '#' | awk '{print $1"\t"$2}' | sort -k1,1 -k2,2n > 166_potato_Allchr.gatk.filter.${i}.pop13.xls

echo -e "american\t${i} private"
./retain_samples_in_vcf.py 166_potato_Allchr.gatk.filter.${i}.vcf pop12.xls 4|grep -v '#' | awk '{print $1"\t"$2}' | sort -k1,1 -k2,2n > 166_potato_Allchr.gatk.filter.${i}.pop12.xls
done


### private variants (redo)
for i in SNP indel
do
echo -e "starch\t${i} private"
sort 166_potato_Allchr.gatk.filter.${i}.pop1_starch.xls -o 166_potato_Allchr.gatk.filter.${i}.pop1_starch.xls
sort 166_potato_Allchr.gatk.filter.${i}.pop23.xls -o 166_potato_Allchr.gatk.filter.${i}.pop23.xls

comm -23 166_potato_Allchr.gatk.filter.${i}.pop1_starch.xls 166_potato_Allchr.gatk.filter.${i}.pop23.xls > 166_potato_Allchr.gatk.filter.${i}.starch.private.xls
wc -l 166_potato_Allchr.gatk.filter.${i}.starch.private.xls

sort 166_potato_Allchr.gatk.filter.${i}.pop2_european.xls -o 166_potato_Allchr.gatk.filter.${i}.pop2_european.xls
sort 166_potato_Allchr.gatk.filter.${i}.pop13.xls -o 166_potato_Allchr.gatk.filter.${i}.pop13.xls

comm -23 166_potato_Allchr.gatk.filter.${i}.pop2_european.xls 166_potato_Allchr.gatk.filter.${i}.pop13.xls > 166_potato_Allchr.gatk.filter.${i}.european.private.xls
wc -l 166_potato_Allchr.gatk.filter.${i}.european.private.xls

sort 166_potato_Allchr.gatk.filter.${i}.pop3_american.xls -o 166_potato_Allchr.gatk.filter.${i}.pop3_american.xls
sort 166_potato_Allchr.gatk.filter.${i}.pop12.xls -o 166_potato_Allchr.gatk.filter.${i}.pop12.xls

comm -23 166_potato_Allchr.gatk.filter.${i}.pop3_american.xls 166_potato_Allchr.gatk.filter.${i}.pop12.xls > 166_potato_Allchr.gatk.filter.${i}.american.private.xls
wc -l 166_potato_Allchr.gatk.filter.${i}.american.private.xls
done

### group specific variants impacting genes (MODERATE/HIGH)
for n in starch european american
do
   for i in SNP indel
   do
	   sort -k1,1 -k2,2n 166_potato_Allchr.gatk.filter.${i}.${n}.private.xls -o 166_potato_Allchr.gatk.filter.${i}.${n}.private.xls
       zcat 166_potato_Allchr.gatk.filter.${i}.vcf.snpeff.gz | ./extract_pos_from_vcf.py 166_potato_Allchr.gatk.filter.${i}.${n}.private.xls > 166_potato_Allchr.gatk.filter.${i}.${n}.private.vcf.snpeff
   done
   cat 166_potato_Allchr.gatk.filter.*.${n}.private.vcf.snpeff | egrep 'MODERATE|HIGH' > 166_potato_Allchr.gatk.filter.${n}.private.vcf.snpeff_MODERATE_HIGH.xls
   cat 166_potato_Allchr.gatk.filter.${n}.private.vcf.snpeff_MODERATE_HIGH.xls|awk '{print $8}'|awk -F 'ANN=' '{print $NF}'|sed 's/,/\n/g'|egrep 'MODERATE|HIGH'|awk -F '|' '{print $7}'|sed 's/-/\n/g'|sort |uniq > 166_potato_Allchr.gatk.filter.${n}.private.vcf.snpeff.gene.xls 
   echo -e "$n\tMODERATE_HIGH_number"
   wc -l 166_potato_Allchr.gatk.filter.${n}.private.vcf.snpeff_MODERATE_HIGH.xls
   echo -e "$n\timpact genes"
   wc -l 166_potato_Allchr.gatk.filter.${n}.private.vcf.snpeff.gene.xls
done

### private gene GO enrichement
./output_go_ipr_from_annot.py DM_v6.1_GO_IPR_Swiss_Ara_KEGG.txt DM_v6_GO.xls2 DM_v6_IPR.xls > DM_v6_GO.xls

grep -v '^$' DM_v6_GO.xls2 | awk 'NF!=1' > t && mv t DM_v6_GO.xls2
grep -v '^$' DM_v6_GO.xls | awk 'NF!=1' > t && mv t DM_v6_GO.xls


for n in starch
do
./go.R DM_v6_GO.xls2 166_potato_Allchr.gatk.filter.${n}.private.vcf.snpeff.gene.xls 166_potato_Allchr.gatk.filter.${n}.private.vcf.snpeff.gene_GO.xls
done


### private geen IPR enrichment

./output_ipr_from_gff3.py DM_v6.1_GO_IPR_Swiss_Ara_KEGG.txt no | awk 'NF!=1' > DM_v6_IPR_gene_infor.xls

awk -F '\t' '{print $6}' <(cat DM_v6.1_GO_IPR_Swiss_Ara_KEGG.txt) |awk '$1!="-" && $1!="IPR"'|sed 's/;/\n/g'|grep -v '^$' | sort | uniq > DM_v6_IPR_infor.xls

rm starch_private_variants_MODERATE_HIGH_gene_IPR.xls
cat 166_potato_Allchr.gatk.filter.starch.private.vcf.snpeff.gene.xls | while read i
do
    grep $i DM_v6_IPR_gene_infor.xls | awk 'NF!=1' >> starch_private_variants_MODERATE_HIGH_gene_IPR.xls
done

./output_ipr_diff.py starch_private_variants_MODERATE_HIGH_gene_IPR.xls DM_v6_IPR_gene_infor.xls DM_v6_IPR_infor.xls > starch_private_variants_MODERATE_HIGH_gene_IPR_ratio_infor.xls


rm tmp fisher_test_p_value.xls
grep -v '#' starch_private_variants_MODERATE_HIGH_gene_IPR_ratio_infor.xls | while read i 
    do
    echo $i | awk '{print $1"\t"$(NF-1)-$1}' >> tmp
    echo $i | awk '{print $3"\t"$(NF)-$3}' >> tmp
    ./fisher_test.R tmp | grep 'p-value' | awk '{print $NF}' >> fisher_test_p_value.xls
    rm tmp
    done
sed -i '1i fisher_p_value' fisher_test_p_value.xls
paste starch_private_variants_MODERATE_HIGH_gene_IPR_ratio_infor.xls fisher_test_p_value.xls > tmp ; mv tmp starch_private_variants_MODERATE_HIGH_gene_IPR_ratio_infor.xls
awk '{print $NF}' starch_private_variants_MODERATE_HIGH_gene_IPR_ratio_infor.xls > tmp
./fdr.R tmp fdr.xls
cat <(echo -e 'fdr_value') fdr.xls > fdr_column.xls
paste starch_private_variants_MODERATE_HIGH_gene_IPR_ratio_infor.xls fdr_column.xls > tmp && sort -t $'\t' -k9,9g tmp > starch_private_variants_MODERATE_HIGH_gene_IPR_ratio_infor.xls

