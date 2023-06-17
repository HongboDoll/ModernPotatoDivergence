#!/bin/bash

ls vernei/ | awk -F '.' '{print $1}' | while read i
do

cat <(echo -e "Chr\tPos\tPos2\tD\tf_d\tf_dM\td_f") <(awk 'NR!=1' vernei/${i}.txt | awk 'BEGIN{OFS="\t"}{if($1<10){print "chr0"$1,$2,$3,$4,$5,$6,$7}else{print "chr"$1,$2,$3,$4,$5,$6,$7}}') > ${i}.xls

cat <(cut -f '1,2,6' ${i}.xls | head -1) <(cut -f '1,2,6' ${i}.xls | awk 'NR!=1' | awk 'BEGIN{OFS="\t"}{if($3<0){print $1,$2,0}else{print $1,$2,$3}}') > t
./plot_introgression.R t starch_introgression_${i}_fdm.pdf 12 chr

done

######## Plot introgression between Altus and Kardent
echo -e "Chr\tStart\tEnd" > DM_v6.1_karyotype.xls
out_len.py DM_v6.1.fa |grep chr | sort -k1,1 | sed 's/>//g' | awk 'BEGIN{OFS="\t"}{print $1,0,$2}' >> DM_v6.1_karyotype.xls

echo -e "Type\tShape\tChr\tStart\tEnd\tcolor" > DM_v6.1_NLR.xls
awk '$3=="mRNA"' Solanum_tuberosumDM.integratedNLR.gff3  |grep NLR |grep -v scaffold|awk 'BEGIN{OFS="\t"}{print "NLR","circle",$1,$4,$5,"8df9a2"}' >> DM_v6.1_NLR.xls

echo -e "Chr\tStart\tEnd\tValue" > DM_v6.1_Altus_vernei_introgression.xls
awk 'NR!=1' European_Altus_vernei_localFstats__500_250.xls | awk 'BEGIN{OFS="\t"}{if($6<0){print $1,$2,$3,0}else{print $1,$2,$3,$6}}' | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4}' >> DM_v6.1_Altus_vernei_introgression.xls

echo -e "Chr\tStart\tEnd\tValue" > DM_v6.1_Kardent_vernei_introgression.xls
awk 'NR!=1' European_Kardent_vernei_localFstats__500_250.xls | awk 'BEGIN{OFS="\t"}{if($6<0){print $1,$2,$3,0}else{print $1,$2,$3,$6}}' | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4}' >> DM_v6.1_Kardent_vernei_introgression.xls


######## compare introgression proportions bewteen yield-related qtls and R genes
cat Solanum_tuberosumDM.integratedNLR.gff3|grep NLR|awk '$3=="mRNA"'|cut -f '1,4,5'|grep -v scaffold|sort -k1,1 -k2,2n |bedtools merge -d 1000000 > Solanum_tuberosumDM.integratedNLR_merge_1Mb.bed

ls *__500_250.xls | grep -v '_Starch_' | while read i
do
awk 'NR!=1' $i|cut -f '1-3,6'|awk 'BEGIN{OFS="\t"}{if($4<0){print $1,$2,$3,0}else{print $1,$2,$3,$4}}' | sort -k1,1 -k2,2n > ${i}_fdM.bed

name=`echo $i | awk -F '_' '{print $2}'`
all=`cat ${i}_fdM.bed | awk '{i+=(($3-$2+1)*$4)}END{print i}'`
# R gene introgression proportion; 732 NLR; 距离1mb以内合并，44614704 bp
r=`bedtools intersect -wa -a ${i}_fdM.bed -b Solanum_tuberosumDM.integratedNLR_merge_1Mb.bed | awk '{i+=(($3-$2+1)*$4)}END{print i/"'"$all"'"*100}'`

# yield-related qtl introgression proportion; 25qtl, tsc/ty/tsy, 169086141 bp
yield=`bedtools intersect -wa -a ${i}_fdM.bed -b reported_qtl_gwas_interval_yield.bed | awk '{i+=(($3-$2+1)*$4)}END{print i/"'"$all"'"*100}'`

# sga-related qtl introgression proportion; 5qtl, sga, 87471277 bp
sga=`bedtools intersect -wa -a ${i}_fdM.bed -b reported_qtl_gwas_interval_sga.bed | awk '{i+=(($3-$2+1)*$4)}END{print i/"'"$all"'"*100}'`
echo -e "${name}\t${sga}\t${r}\t${yield}"
done

### all 23 starch 
awk 'NR!=1' European_Starch_vernei_localFstats__500_250.xls|cut -f '1-3,6'|awk 'BEGIN{OFS="\t"}{if($4<0){print $1,$2,$3,0}else{print $1,$2,$3,$4}}' | sort -k1,1 -k2,2n > European_Starch_vernei_localFstats__500_250.xls_fdM.bed
i=European_Starch_vernei_localFstats__500_250.xls
all=`cat ${i}_fdM.bed | awk '{i+=(($3-$2+1)*$4)}END{print i}'`
### 6.35 Mb of the 44.6 Mb NLR regions
r=`bedtools intersect -wa -a ${i}_fdM.bed -b Solanum_tuberosumDM.integratedNLR_merge_1Mb.bed | awk '{i+=(($3-$2+1)*$4)}END{print i/"'"$all"'"*100}'`

### 8.56 Mb of the 169.1 MB qtl regions
# yield-related qtl introgression proportion; 25qtl, tsc/ty/tsy, 169086141 bp
yield=`bedtools intersect -wa -a ${i}_fdM.bed -b reported_qtl_gwas_interval_yield.bed | awk '{i+=(($3-$2+1)*$4)}END{print i/"'"$all"'"*100}'`
echo -e "23_starch\t${r}\t${yield}"

### 3.81 Mb of the 87.5 Mb qtl regions
# sga-related qtl introgression proportion; 5qtl, sga, 87471277 bp
sga=`bedtools intersect -wa -a ${i}_fdM.bed -b reported_qtl_gwas_interval_sga.bed | awk '{i+=(($3-$2+1)*$4)}END{print i/"'"$all"'"*100}'`
echo -e "23_starch\t${sga}\t${r}\t${yield}"

