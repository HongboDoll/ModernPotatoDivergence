#!/bin/bash

cd /media/scratchpad_02/li322/work/04_137_atlas_population_analyses/nucleotide_diversity_fst/qtl_gwas_gene_overlap

fst_cutoff=0.2093
pai_cutoff=1.84654

#cut -f '1-5' reported_qtl_gwas_interval.xls  |sort -k1,1 -k2,2n > reported_qtl_gwas_interval.bed
#
#awk 'NR!=1' 166_potato_european_vs_starch_nucleotide_diversity_region.xls | awk '$4>"'"$pai_cutoff"'"' | sed 's/ /\t/g' | sort -k1,1 -k2,2n | bedtools merge -d 50000 > 166_potato_european_vs_starch_nucleotide_diversity_region.bed
#
#sort 166_potato_starch_vs_european_Fst_region.xls > 166_potato_starch_vs_european_Fst_region.bed
#
#bedtools intersect -a reported_qtl_gwas_interval.bed -b 166_potato_starch_vs_european_Fst_region.bed -u | sort > reported_qtl_gwas_interval_overlapped_fst.xls
#
#bedtools intersect -a reported_qtl_gwas_interval.bed -b 166_potato_european_vs_starch_nucleotide_diversity_region.bed -u | sort > reported_qtl_gwas_interval_overlapped_pai.xls
#
#comm -12 reported_qtl_gwas_interval_overlapped_fst.xls reported_qtl_gwas_interval_overlapped_pai.xls > reported_qtl_gwas_interval_overlapped_both_fst_pai.xls

for n in fst pai
do
awk '{print $1}' reported_qtl_gwas_interval_overlapped_${n}.xls  |sort |uniq > chr.xls
awk 'NR!=1' 166_potato_european_vs_starch_nucleotide_diversity_region.xls | awk '{print $1}' | sort | uniq > chrall.xls
comm -23 chrall.xls chr.xls > chrlack.xls
cat chrlack.xls | while read i
do
	echo -e "${i}\t0\t0\tNA\tNA" >> reported_qtl_gwas_interval_overlapped_${n}.xls
done
sort -k1,1 -k2,2n reported_qtl_gwas_interval_overlapped_${n}.xls -o reported_qtl_gwas_interval_overlapped_${n}.xls
done

cat <(cut -f '1,2,4' 166_potato_3pop_paired_Fst.xls | head -1) <(cut -f '1,2,4' 166_potato_3pop_paired_Fst.xls | awk 'NR!=1' | awk 'BEGIN{OFS="\t"}{if($3<0){print $1,$2,0}else{print $1,$2,$3}}') > t
./plot_fst_qtl.R t 166_potato_starch_vs_european_Fst.pdf 12 chr 0.2093 reported_qtl_gwas_interval_overlapped_fst.xls

cat <(echo -e "Chr\tPos\tPai") <(awk 'NR!=1' 166_potato_3pop_nucleotide_diversity.xls |awk '$5>0.0002&&$4>0.0001'|awk '{if($4!=0){print $1,$2,$5/$4}else{print $1,$2,0}}') > t
./plot_pai_qtl.R t 166_potato_european_vs_starch_nucleotide_diversity.pdf 12 chr 1.84654 reported_qtl_gwas_interval_overlapped_pai.xls

#### local pai for example

cat <(echo -e "Chr\tPos\tPai1\tPai2") <(awk 'NR!=1' 166_potato_3pop_nucleotide_diversity.xls |awk '$5>0.0002&&$4>0.0001'|awk '{print $1,$2,$4,$5}')|sed -n '4983,4993p' > chr05_gwd_pai_4_plot.xls
./plot_local_pai.R chr05_gwd_pai_4_plot.xls chr05_gwd_pai_4_plot.pdf

cat <(echo -e "Chr\tPos\tPai1\tPai2") <(awk 'NR!=1' 166_potato_3pop_nucleotide_diversity.xls |awk '$5>0.0002&&$4>0.0001'|awk '{print $1,$2,$4,$5}')|sed -n '10719,10725p' > chr09_susy3_pai_4_plot.xls
./plot_local_pai.R chr09_susy3_pai_4_plot.xls chr09_susy3_pai_4_plot.pdf

cat <(echo -e "Chr\tPos\tPai1\tPai2") <(awk 'NR!=1' 166_potato_3pop_nucleotide_diversity.xls |awk '$5>0.0002&&$4>0.0001'|awk '{print $1,$2,$4,$5}')|sed -n '13449,13459p' > chr12_susy4_pai_4_plot.xls
./plot_local_pai.R chr12_susy4_pai_4_plot.xls chr12_susy4_pai_4_plot.pdf

