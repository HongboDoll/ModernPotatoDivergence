#!/bin/bash

cd /media/scratchpad_02/li322/work/05_introgression

#gzip -c -d 166_potato_Allchr.gatk.filter.SNP.filterMissMAF.diploid_genotype.vcf.gz > 166_potato_Allchr.gatk.filter.SNP.filterMissMAF.diploid_genotype.vcf

vcf=166_potato_Allchr.gatk.filter.SNP.filterMissMAF.diploid_genotype.vcf
retain_given_sample_in_vcf.py $vcf pop12.xls 2 | gzip > 166_potato_23_starch__47_European_diploid.vcf.gz
retain_given_sample_in_vcf.py $vcf pop2_european_4_introgression.xls 2 | gzip > 166_potato_47_European_diploid.vcf.gz

###### Dsuite

head -1000 367_all_chr_snp_indel_keepPG.vcf | grep '#' | grep -v 'scaffold' > head
cat head <(paste <(grep -v '#' 367_all_chr_snp_indel_keepPG.vcf | grep chr | awk 'length($4)==1&&length($5)==1' | cut -f '1-2') <(grep -v '#' 367_all_chr_snp_indel_keepPG.vcf | grep chr | awk 'length($4)==1&&length($5)==1' | awk '{print "SNP367_"$1"_"$2}') <(grep -v '#' 367_all_chr_snp_indel_keepPG.vcf | grep chr | awk 'length($4)==1&&length($5)==1' | cut -f '4-9999')) > 367_all_chr_snp_indel_keepPG_SNP.vcf

cut -f '1' wild_outgroup.xls > retain
retain_given_sample_in_vcf.py 367_all_chr_snp_indel_keepPG_SNP.vcf retain 2 > 367_all_chr_snp_indel_wild_outgroup.vcf

cp /home/huyong/tmp/166_potato_23_starch__47_European_diploid.vcf.gz 166_potato_23_starch_47_European_diploid.vcf.gzz
gzip -c -d 166_potato_23_starch_47_European_diploid.vcf.gzz | bgzip > 166_potato_23_starch_47_European_diploid.vcf.gz

tabix -f -p vcf 166_potato_23_starch_47_European_diploid.vcf.gz

bgzip 367_all_chr_snp_indel_wild_outgroup.vcf
tabix -f -p vcf 367_all_chr_snp_indel_wild_outgroup.vcf.gz

zcat 166_potato_23_starch_47_European_diploid.vcf.gz | head -100| grep '#' > head

cat head2 <(zcat 367_all_chr_snp_indel_wild_outgroup.vcf.gz | grep -v '#') | ./retain_only_GT_vcf.py | bgzip > t && mv t 367_all_chr_snp_indel_wild_outgroup.vcf.gz
tabix -f -p vcf 367_all_chr_snp_indel_wild_outgroup.vcf.gz
cat head <(zcat 166_potato_23_starch_47_European_diploid.vcf.gz | grep -v '#') | ./retain_only_GT_vcf.py | bgzip > t && mv t 166_potato_23_starch_47_European_diploid.vcf.gz
tabix -f -p vcf 166_potato_23_starch_47_European_diploid.vcf.gz

bcftools merge -0  166_potato_23_starch_47_European_diploid.vcf.gz 367_all_chr_snp_indel_wild_outgroup.vcf.gz > 367_snp_wild_outgroup_23_starch_47_European.vcf


vcf=367_snp_wild_outgroup_23_starch_47_European.vcf

plink --vcf $vcf --maf 0.05 --geno 0.2 --recode vcf-iid --allow-no-sex --biallelic-only -out 367_snp_wild_outgroup_23_starch_47_European_4_ABBA --allow-extra-chr

cat <(cat wild_outgroup.xls|grep Outgroup | cut -f '1,3') <(grep -v Outgroup wild_outgroup.xls | cut -f '1,2')|sort -k2,2 > wild_outgroup_set.xls

cat Starch_European_set.xls  wild_outgroup_set.xls > set.txt

out_len.py DM_v6.1.fa | grep chr | sed 's/>//g' | sed 's/chr0//g' | sed 's/chr//g' > DM_v6.1_len.fa
bedtools makewindows -g DM_v6.1_len.fa -w 100000 > DM_v6.1_window_100k.xls

head -100 367_snp_wild_outgroup_23_starch_47_European_4_ABBA.vcf | grep '#' > head
rm -rf Dtrios_results; mkdir Dtrios_results
cd Dtrios_results
ln -s ../367_snp_wild_outgroup_23_starch_47_European_4_ABBA.vcf; ln -s ../head; ln -s ../set.txt
for i in `seq 1 12`
do
	echo """#!/bin/bash
	cat head <(awk '\$1==$i' 367_snp_wild_outgroup_23_starch_47_European_4_ABBA.vcf) > 367_snp_wild_outgroup_23_starch_47_European_4_ABBA_chr${i}.vcf
	Dsuite Dtrios  -o output_chr${i} 367_snp_wild_outgroup_23_starch_47_European_4_ABBA_chr${i}.vcf set.txt
	""" > chr${i}.sh && chmod 755 chr${i}.sh
	sbatch -J chr${i} -p queue1 --qos=queue1 -N 1 --ntasks-per-node=1 -e %x.err -o %x.out "./chr${i}.sh"
done
cd -

Dsuite DtriosCombine -o whole_genome Dtrios_results/output_chr1 Dtrios_results/output_chr2 Dtrios_results/output_chr3 Dtrios_results/output_chr4 Dtrios_results/output_chr5 Dtrios_results/output_chr6 Dtrios_results/output_chr7 Dtrios_results/output_chr8 Dtrios_results/output_chr9 Dtrios_results/output_chr10 Dtrios_results/output_chr11 Dtrios_results/output_chr12


rm -rf Dinvestigate_results; mkdir Dinvestigate_results
cd Dinvestigate_results
ln -s ../set.txt; ln -s ../367_snp_wild_outgroup_23_starch_47_European_4_ABBA.vcf
grep PG set.txt |grep -v Outgroup | awk '{print $2}'|sort | uniq | while read sp
do
    echo -e "European\tStarch\t${sp}" > test_trios_${sp}.txt
	echo -e """#!/bin/bash
	Dsuite  Dinvestigate -w 500,250 367_snp_wild_outgroup_23_starch_47_European_4_ABBA.vcf set.txt test_trios_${sp}.txt
	""" > ${sp}.sh && chmod 755 ${sp}.sh
	sbatch -J $sp -p gpu-3090 -N 1 --ntasks-per-node=1 -e %x.err -o %x.out "./${sp}.sh" ### --nodelist 指定节点运行  login7 login8好像有问题，至少跑NLR-tracker不行
done
cd -


######### statistics
./format_order_whole_genome_D_Z.py whole_genome_combined_BBAA.txt wild_order.xls > whole_genome_combined_D_Z_order.xls

#rm per_chromosome_D.xls
echo -e "chr\t\c" >> per_chromosome_D.xls
awk -F '\t' '{print $1}' whole_genome_combined_D_Z_order.xls | while read i
do
	echo -e "${i}\t\c" >> per_chromosome_D.xls
done
echo -e "" >> per_chromosome_D.xls

for i in `seq 1 12`
do
	./format_per_chromosome_D.py Dtrios/output_chr${i}_BBAA.txt whole_genome_combined_D_Z_order.xls chr${i} >> per_chromosome_D.xls
done

###### overall introgressed ratio
awk 'NR!=1' European_Starch_vernei_localFstats__500_250.xls | awk 'BEGIN{OFS="\t"}{if($6<0){print $1,$2,$3,0}else{print $1,$2,$3,$6}}' | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4}'|awk '{i+=(($3-$2)*$4)}END{print i}'

###### plot introgression regions

ls Dinvestigate_500_250/ | awk -F '.' '{print $1}' | while read i
do

cat <(echo -e "Chr\tPos\tPos2\tD\tf_d\tf_dM\td_f") <(awk 'NR!=1' Dinvestigate_500_250/${i}.txt | awk 'BEGIN{OFS="\t"}{if($1<10){print "chr0"$1,$2,$3,$4,$5,$6,$7}else{print "chr"$1,$2,$3,$4,$5,$6,$7}}') > ${i}.xls

cat <(cut -f '1,2,6' ${i}.xls | head -1) <(cut -f '1,2,6' ${i}.xls | awk 'NR!=1' | awk 'BEGIN{OFS="\t"}{if($3<0){print $1,$2,0}else{print $1,$2,$3}}') > t
./plot_introgression.R t starch_introgression_${i}_fdm.pdf 12 chr

done

###### Plot with yield qtls

awk '{print $1}' reported_qtl_gwas_interval_yield.bed  |sort |uniq > chr.xls
awk 'NR!=1' t | awk '{print $1}' | sort | uniq > chrall.xls
comm -23 chrall.xls chr.xls > chrlack.xls
cat chrlack.xls | while read i
do
	echo -e "${i}\t0\t0\tNA\tNA" >> reported_qtl_gwas_interval_yield.bed
done
sort -k1,1 -k2,2n reported_qtl_gwas_interval_yield.bed -o reported_qtl_gwas_interval_yield.bed

i=European_Starch_vernei_localFstats__500_250

cat <(cut -f '1,2,6' ${i}.xls | head -1) <(cut -f '1,2,6' ${i}.xls | awk 'NR!=1' | awk 'BEGIN{OFS="\t"}{if($3<0){print $1,$2,0}else{print $1,$2,$3}}') > t
./plot_introgression_qtl.R t starch_introgression_${i}_fdm_QTL.pdf 12 chr reported_qtl_gwas_interval_yield.bed

###### overlap with qtl
introgression_cutoff=0.177504

awk 'NR!=1' European_Starch_vernei_localFstats__50_25.xls | awk '$6>"'"$introgression_cutoff"'"' | sed 's/ /\t/g' | sort -k1,1 -k2,2n | bedtools merge -d 50000 > European_Starch_vernei_localFstats__50_25.bed

bedtools intersect -a reported_qtl_gwas_interval.bed -b European_Starch_vernei_localFstats__50_25.bed -u | sort > reported_qtl_gwas_interval_overlapped_introgression.xls

for n in introgression
do
awk '{print $1}' reported_qtl_gwas_interval_overlapped_${n}.xls  |sort |uniq > chr.xls
awk 'NR!=1' European_Starch_vernei_localFstats__50_25.xls | awk '{print $1}' | sort | uniq > chrall.xls
comm -23 chrall.xls chr.xls > chrlack.xls
cat chrlack.xls | while read i
do
	echo -e "${i}\t0\t0\tNA\tNA" >> reported_qtl_gwas_interval_overlapped_${n}.xls
done
sort -k1,1 -k2,2n reported_qtl_gwas_interval_overlapped_${n}.xls -o reported_qtl_gwas_interval_overlapped_${n}.xls
done


