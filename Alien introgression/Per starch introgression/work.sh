#!/bin/bash

####### Dsuite
awk '{if($2=="Starch"){print $1"\t"$1}else{print $0}}' set.txt > set.txt_per_starch
rm -rf Dinvestigate_results; mkdir Dinvestigate_results
cd Dinvestigate_results
ln -s ../set.txt_per_starch; ln -s ../367_snp_wild_outgroup_23_starch_47_European_4_ABBA.vcf
grep PG set.txt_per_starch |grep -v Outgroup | awk '{print $2}'|sort | uniq | while read sp
do
	grep -v PG set.txt_per_starch|grep -v European | awk '{print $2}'|sort | uniq | while read starch
	do
    echo -e "European\t${starch}\t${sp}" > test_trios_${starch}_${sp}.txt
	echo -e """#!/bin/bash
	Dsuite  Dinvestigate -w 500,250 367_snp_wild_outgroup_23_starch_47_European_4_ABBA.vcf set.txt_per_starch test_trios_${starch}_${sp}.txt
	""" > ${starch}_${sp}.sh && chmod 755 ${starch}_${sp}.sh
	sbatch -J ${starch}_$sp -p gpu-3090 -N 1 --ntasks-per-node=1 -e %x.err -o %x.out "./${starch}_${sp}.sh" ### --nodelist 指定节点运行  login7 login8好像有问题，至少跑NLR-tracker不行
	done
done
cd -

mkdir vernei
cp Dinvestigate_results/*vernei*localFstats*txt vernei

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

ls *__500_250.xls | while read i
do
awk 'NR!=1' $i|cut -f '1-3,6'|awk 'BEGIN{OFS="\t"}{if($4<0){print $1,$2,$3,0}else{print $1,$2,$3,$4}}' | sort -k1,1 -k2,2n > ${i}_fdM.bed

name=`echo $i | awk -F '_' '{print $2}'`
all=`cat ${i}_fdM.bed | awk '{i+=(($3-$2+1)*$4)}END{print i}'`
# R gene introgression proportion; 732 NLR; 距离1mb以内合并，44614704 bp
r=`bedtools intersect -wa -a ${i}_fdM.bed -b Solanum_tuberosumDM.integratedNLR_merge_1Mb.bed | awk '{i+=(($3-$2+1)*$4)}END{print i/"'"$all"'"*100}'`

# yield-related qtl introgression proportion; 25qtl, tsc/ty/tsy, 169086141 bp
yield=`bedtools intersect -wa -a ${i}_fdM.bed -b reported_qtl_gwas_interval_yield.bed | awk '{i+=(($3-$2+1)*$4)}END{print i/"'"$all"'"*100}'`

echo -e "${name}\t${r}\t${yield}"
done

### all 23 starch 
awk 'NR!=1' European_Starch_vernei_localFstats__500_250.xls|cut -f '1-3,6'|awk 'BEGIN{OFS="\t"}{if($4<0){print $1,$2,$3,0}else{print $1,$2,$3,$4}}' | sort -k1,1 -k2,2n > European_Starch_vernei_localFstats__500_250.xls_fdM.bed
i=European_Starch_vernei_localFstats__500_250.xls
all=`cat ${i}_fdM.bed | awk '{i+=(($3-$2+1)*$4)}END{print i}'`
r=`bedtools intersect -wa -a ${i}_fdM.bed -b Solanum_tuberosumDM.integratedNLR_merge_1Mb.bed | awk '{i+=(($3-$2+1)*$4)}END{print i/"'"$all"'"*100}'`

# yield-related qtl introgression proportion; 25qtl, tsc/ty/tsy, 169086141 bp
yield=`bedtools intersect -wa -a ${i}_fdM.bed -b reported_qtl_gwas_interval_yield.bed | awk '{i+=(($3-$2+1)*$4)}END{print i/"'"$all"'"*100}'`
echo -e "23_starch\t${r}\t${yield}"

