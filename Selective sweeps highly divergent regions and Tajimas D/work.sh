#!/bin/bash

cd /media/scratchpad_02/li322/work/04_137_atlas_population_analyses/nucleotide_diversity_fst

######## 67005864 SNP with MAF 0.05 and missing call rate 0.2 filtered are used 
######## 8375 SNP per 100kb
#for i in {01..12}
#do
#    echo -e """#!/bin/bash
#    cd /media/scratchpad_02/li322/work/04_137_atlas_population_analyses/nucleotide_diversity_fst
#	cd filter_vcf/chr${i}
#	bgzip 166_potato_Allchr.gatk.filter.SNP.filterMissMAF.chr${i}.vcf
#	tabix -f -p vcf 166_potato_Allchr.gatk.filter.SNP.filterMissMAF.chr${i}.vcf.gz
#	cd -
#    """ > vcf_chr${i}.sh && chmod 755 vcf_chr${i}.sh
#    qsub -cwd -N chr${i} -q stat.q -j y -V -b y -S /bin/bash "//media/scratchpad_02/li322/work/04_137_atlas_population_analyses/nucleotide_diversity_fst/vcf_chr${i}.sh"
#done

window=100000   # argv[5]
slide=50000   # argv[5]/2, half of the window size
gff=DM_v6.1.gff3   # argv[4]

### To read in SNP data from a subset of individuals the parameter samplenames requires an character vector including the indi
### vcf_handle  <- .Call("VCF_open",filename)
### ind         <- .Call("VCF_getSampleNames",vcf_handle)
### samplenames <- ind[1:10]
### pop1_starch.txt: starch variety
### pop2_european.xls: EU variety
### pop3_american.xls: American variety (US and CA)

rm -rf 166_potato_pop_ana_results; mkdir 166_potato_pop_ana_results
grep '>chr' DM_v6.1.fa | sed 's/>//g' | while read i
do
       mkdir -p 166_potato_pop_ana_results/${i}
       end=`out_len.py DM_v6.1.fa | grep $i | awk '{print $2}'`  # argv[3]
       ref=$i   # argv[2]
       nuc_div2=${i}_nucleotide_diversity.xls2   # argv[6]
       fst2=${i}_Fst.xls2   # argv[7]
       tajimasd2=${i}_TajimasD.xls2   # argv[8]

       nuc_div=${i}_nucleotide_diversity.xls   # argv[6]
       fst=${i}_Fst.xls   # argv[7]
       tajimasd=${i}_TajimasD.xls   # argv[8]
       vcf=$PWD/filter_vcf/${i}/166_potato_Allchr.gatk.filter.SNP.filterMissMAF.${i}.vcf.gz   # argv[1]

	   fst3=${i}_Fst_each_site.xls3
       tajimasd3=${i}_TajimasD_diploid.xls3
	   echo -e """#!/bin/bash
	   
	   cd /media/scratchpad_02/li322/work/04_137_atlas_population_analyses/nucleotide_diversity_fst/166_potato_pop_ana_results/${i}/
	   
       ../../popgenome.R $vcf $ref $end ../../$gff $window $nuc_div $fst $tajimasd
       ../../popgenome.R2 $vcf $ref $end ../../$gff $window $nuc_div2 $fst2 $tajimasd2  #### one vs. the othres
       ../../popgenome.R3 $vcf $ref $end ../../$gff 1 $fst3 #### FST per site
       ../../popgenome.R4 $vcf $ref $end ../../$gff $window $tajimasd3
	   """ > ${i}_pop_ana.sh && chmod 755 ${i}_pop_ana.sh
		qsub -cwd -N ${i}_pai_fst -pe sharedmem 1 -q seq.q -j y -V -b y -S /bin/bash "//media/scratchpad_02/li322/work/04_137_atlas_population_analyses/nucleotide_diversity_fst/${i}_pop_ana.sh"
done

####### summary

window=100000
slide=50000

grep '>chr' DM_v6.1.fa | sed 's/>//g' | while read i
do

       ref=$i   # argv[2]
       nuc_div=166_potato_pop_ana_results/${i}/${i}_nucleotide_diversity.xls   # argv[6]
       fst=166_potato_pop_ana_results/${i}/${i}_Fst.xls   # argv[7]
       tajimasd=166_potato_pop_ana_results/${i}/${i}_TajimasD.xls   # argv[8]
	   
       nuc_div2=166_potato_pop_ana_results/${i}/${i}_nucleotide_diversity.xls2   # argv[6]
       fst2=166_potato_pop_ana_results/${i}/${i}_Fst.xls2   # argv[7]
       tajimasd2=166_potato_pop_ana_results/${i}/${i}_TajimasD.xls2   # argv[8]
	   
	   awk 'NR!=1' 166_potato_pop_ana_results/${i}/${i}_nucleotide_diversity.xls | awk -F '\t' 'BEGIN{OFS="\t"}{print "'"$i"'",($1-1)*"'"$slide"'",($1-1)*"'"$slide"'"+"'"$window"'",$2,$3,$4}' > t && mv t 166_potato_pop_ana_results/${i}/${i}_nucleotide_diversity.xls
	   awk 'NR!=1' 166_potato_pop_ana_results/${i}/${i}_nucleotide_diversity.xls2 | awk -F '\t' 'BEGIN{OFS="\t"}{print "'"$i"'",($1-1)*"'"$slide"'",($1-1)*"'"$slide"'"+"'"$window"'",$2,$3}' > t && mv t 166_potato_pop_ana_results/${i}/${i}_nucleotide_diversity.xls2
	   awk 'NR!=1' 166_potato_pop_ana_results/${i}/${i}_TajimasD.xls | awk -F '\t' 'BEGIN{OFS="\t"}{print "'"$i"'",($1-1)*"'"$slide"'",($1-1)*"'"$slide"'"+"'"$window"'",$2,$3,$4}' > t && mv t 166_potato_pop_ana_results/${i}/${i}_TajimasD.xls
	   awk 'NR!=1' 166_potato_pop_ana_results/${i}/${i}_Fst.xls2 | awk -F '\t' 'BEGIN{OFS="\t"}{print "'"$i"'",($1-1)*"'"$slide"'",($1-1)*"'"$slide"'"+"'"$window"'",$2}' > t && mv t 166_potato_pop_ana_results/${i}/${i}_Fst.xls2
	   awk 'NR!=1' 166_potato_pop_ana_results/${i}/${i}_Fst.xls | awk -F '\t' 'BEGIN{OFS="\t"}{print "'"$i"'",($1-1)*"'"$slide"'",($1-1)*"'"$slide"'"+"'"$window"'",$2,$3,$4}' > t && mv t 166_potato_pop_ana_results/${i}/${i}_Fst.xls

done

cat <(echo -e "Chr\tPos\tPos2\tPop1_starch\tPop2_european\tPop3_american") <(cat 166_potato_pop_ana_results/*/*_nucleotide_diversity.xls) > 166_potato_3pop_nucleotide_diversity.xls
cat <(echo -e "Chr\tPos\tPos2\tPop1_starch\tPop23_european_american") <(cat 166_potato_pop_ana_results/*/*_nucleotide_diversity.xls2) > 166_potato_starch_vs_others_nucleotide_diversity.xls
cat <(echo -e "Chr\tPos\tPos2\tPop1_starch\tPop2_european\tPop3_american") <(cat 166_potato_pop_ana_results/*/*_TajimasD.xls) > 166_potato_3pop_TajimasD.xls
cat <(echo -e "Chr\tPos\tPos2\tStarch_vs_european\tStarch_vs_american\tEuropean_vs_american") <(cat 166_potato_pop_ana_results/*/*_Fst.xls) > 166_potato_3pop_paired_Fst.xls
cat <(echo -e "Chr\tPos\tPos2\tPop1_starch_vs_Pop23_european_american") <(cat 166_potato_pop_ana_results/*/*_Fst.xls2) > 166_potato_starch_vs_others_Fst.xls


####### Plot

####### FST
cat <(cut -f '1,2,4' 166_potato_starch_vs_others_Fst.xls | head -1) <(cut -f '1,2,4' 166_potato_starch_vs_others_Fst.xls | awk 'NR!=1' | awk 'BEGIN{OFS="\t"}{if($3<0){print $1,$2,0}else{print $1,$2,$3}}') > t
./plot_fst.R t 166_potato_starch_vs_others_Fst.pdf 12 chr 0.2135

cat <(cut -f '1,2,4' 166_potato_3pop_paired_Fst.xls | head -1) <(cut -f '1,2,4' 166_potato_3pop_paired_Fst.xls | awk 'NR!=1' | awk 'BEGIN{OFS="\t"}{if($3<0){print $1,$2,0}else{print $1,$2,$3}}') > t
./plot_fst.R t 166_potato_starch_vs_european_Fst.pdf 12 chr 0.2093

cat <(cut -f '1,2,5' 166_potato_3pop_paired_Fst.xls | head -1) <(cut -f '1,2,5' 166_potato_3pop_paired_Fst.xls | awk 'NR!=1' | awk 'BEGIN{OFS="\t"}{if($3<0){print $1,$2,0}else{print $1,$2,$3}}') > t
./plot_fst.R t 166_potato_starch_vs_american_Fst.pdf 12 chr 0.2760

cat <(cut -f '1,2,6' 166_potato_3pop_paired_Fst.xls | head -1) <(cut -f '1,2,6' 166_potato_3pop_paired_Fst.xls | awk 'NR!=1' | awk 'BEGIN{OFS="\t"}{if($3<0){print $1,$2,0}else{print $1,$2,$3}}') > t
./plot_fst.R t 166_potato_european_vs_american_Fst.pdf 12 chr 0.1323

###### Tajima's D
cat <(cut -f '1,2,4' 166_potato_3pop_TajimasD.xls | head -1) <(cut -f '1,2,4' 166_potato_3pop_TajimasD.xls | awk 'NR!=1' | sed 's/NA/0/g') > t
./plot_tajimasD.R t 166_potato_starch_TajimasD.pdf 12 chr

cat <(cut -f '1,2,5' 166_potato_3pop_TajimasD.xls | head -1) <(cut -f '1,2,5' 166_potato_3pop_TajimasD.xls | awk 'NR!=1' | sed 's/NA/0/g') > t
./plot_tajimasD.R t 166_potato_european_TajimasD.pdf 12 chr

####### Pai
cat <(echo -e "Chr\tPos\tPai") <(awk 'NR!=1' 166_potato_3pop_nucleotide_diversity.xls |awk '$5>0.0002&&$4>0.0001'|awk '{if($5!=0){print $1,$2,$4/$5}else{print $1,$2,0}}') > t
./plot_pai.R t 166_potato_starch_vs_european_nucleotide_diversity.pdf 12 chr 2.62225

cat <(echo -e "Chr\tPos\tPai") <(awk 'NR!=1' 166_potato_3pop_nucleotide_diversity.xls |awk '$5>0.0002&&$4>0.0001'|awk '{if($4!=0){print $1,$2,$5/$4}else{print $1,$2,0}}') > t
cat <(echo -e "Chr\tPos\tPai") <(awk 'NR!=1' 166_potato_3pop_nucleotide_diversity.xls |awk '$5>0.0002&&$4>0.0001'|awk '{if($4!=0){print $1,$2,$3,$5/$4}else{print $1,$2,$3,0}}') > 166_potato_european_vs_starch_nucleotide_diversity_region.xls

./plot_pai.R t 166_potato_european_vs_starch_nucleotide_diversity.pdf 12 chr 1.84654


#cat <(echo -e "Chr\tPos\tPai") <(awk 'NR!=1' 166_potato_starch_vs_others_nucleotide_diversity.xls |awk '{if($5!=0){print $1,$2,$4/$5}else{print $1,$2,0}}') > t
#./plot_pai_log10.R t 166_potato_starch_vs_others_nucleotide_diversity.pdf 12 chr 2.46
#
#cat <(echo -e "Chr\tPos\tPai") <(awk 'NR!=1' 166_potato_starch_vs_others_nucleotide_diversity.xls |awk '{if($4!=0){print $1,$2,$5/$4}else{print $1,$2,0}}') > t
#./plot_pai_log10.R t 166_potato_others_vs_starch_nucleotide_diversity.pdf 12 chr 2.8471

cat <(echo -e "Chr\tPos\tPai") <(awk 'NR!=1' 166_potato_3pop_nucleotide_diversity.xls |awk '{if($6!=0){print $1,$2,$4/$6}else{print $1,$2,0}}') > t
./plot_pai_log10.R t 166_potato_starch_vs_american_nucleotide_diversity.pdf 12 chr 3.36

cat <(echo -e "Chr\tPos\tPai") <(awk 'NR!=1' 166_potato_3pop_nucleotide_diversity.xls |awk '{if($6!=0){print $1,$2,$5/$6}else{print $1,$2,0}}') > t
./plot_pai_log10.R t 166_potato_european_vs_american_nucleotide_diversity.pdf 12 chr 1.91

cat <(echo -e "Chr\tPos\tPai") <(awk 'NR!=1' 166_potato_3pop_nucleotide_diversity.xls |awk '{if($4!=0){print $1,$2,$6/$4}else{print $1,$2,0}}') > t
./plot_pai_log10.R t 166_potato_american_vs_starch_nucleotide_diversity.pdf 12 chr 3.68769

cat <(echo -e "Chr\tPos\tPai") <(awk 'NR!=1' 166_potato_3pop_nucleotide_diversity.xls |awk '{if($5!=0){print $1,$2,$6/$5}else{print $1,$2,0}}') > t
./plot_pai_log10.R t 166_potato_american_vs_european_nucleotide_diversity.pdf 12 chr 2.35374

###### local pai plot
sed -n '40,80p' 166_potato_starch_vs_others_nucleotide_diversity.xls | cut -f '1,2,4,5' > t
./plot_local_pai.py t 166_potato_starch_vs_others_nucleotide_diversity_local1.pdf


##### FST stat

cat <(cut -f '1,2,4' 166_potato_3pop_paired_Fst.xls | head -1) <(cut -f '1,2,4' 166_potato_3pop_paired_Fst.xls | awk 'NR!=1' | awk 'BEGIN{OFS="\t"}{if($3<0){print $1,$2,0}else{print $1,$2,$3}}') > 166_potato_starch_vs_european_Fst.xls
#
cat <(cut -f '1,2,5' 166_potato_3pop_paired_Fst.xls | head -1) <(cut -f '1,2,5' 166_potato_3pop_paired_Fst.xls | awk 'NR!=1' | awk 'BEGIN{OFS="\t"}{if($3<0){print $1,$2,0}else{print $1,$2,$3}}') > 166_potato_starch_vs_american_Fst.xls
#
cat <(cut -f '1,2,6' 166_potato_3pop_paired_Fst.xls | head -1) <(cut -f '1,2,6' 166_potato_3pop_paired_Fst.xls | awk 'NR!=1' | awk 'BEGIN{OFS="\t"}{if($3<0){print $1,$2,0}else{print $1,$2,$3}}') > 166_potato_european_vs_american_Fst.xls

### Number of regions
grep -v Pos 166_potato_starch_vs_european_Fst.xls|awk 'BEGIN{OFS="\t"}{if($4>0.2093){print $1,$2,$3}}'|bedtools merge -d 50000 > 166_potato_starch_vs_european_Fst_region.xls

grep -v Pos 166_potato_starch_vs_american_Fst.xls|awk 'BEGIN{OFS="\t"}{if($4>0.2760){print $1,$2,$3}}'|bedtools merge -d 50000 > 166_potato_starch_vs_american_Fst_region.xls

grep -v Pos 166_potato_european_vs_american_Fst.xls|awk 'BEGIN{OFS="\t"}{if($4>0.1323){print $1,$2,$3}}'|bedtools merge -d 50000 > 166_potato_european_vs_american_Fst_region.xls

### length
awk '{i+=($3-$2)}END{print i/1000000}' 166_potato_starch_vs_european_Fst_region.xls
awk '{i+=($3-$2)}END{print i/1000000}' 166_potato_starch_vs_american_Fst_region.xls
awk '{i+=($3-$2)}END{print i/1000000}' 166_potato_european_vs_american_Fst_region.xls


## overlapping genes FST

awk 'NR!=1' DM_v6.1_GO_IPR_Swiss_Ara_KEGG.txt |cut -f '1-3' | sort -k1,1 -k2,2n > DM_v6.1_gene_sort.xls

for i in starch_vs_european starch_vs_american european_vs_american
do
	sort -k1,1 -k2,2n 166_potato_${i}_Fst_region.xls -o 166_potato_${i}_Fst_region.xls
	bedtools intersect -a 166_potato_${i}_Fst_region.xls -b DM_v6.1_gene_sort.xls -wb | awk '{print $4,$5,$6}' > 166_potato_${i}_Fst_region_overlap_gene.xls
	rm 166_potato_${i}_Fst_region_overlap_gene_information.xls
	cat 166_potato_${i}_Fst_region_overlap_gene.xls | while read l
	do
		chr=`echo $l | awk '{print $1}'`
		start=`echo $l | awk '{print $2}'`
		end=`echo $l | awk '{print $3}'`
		grep $chr DM_v6.1_GO_IPR_Swiss_Ara_KEGG.txt | grep $start | grep $end >> 166_potato_${i}_Fst_region_overlap_gene_information.xls
	done
done

## overlapping genes pai
awk '$4>1.8465' 166_potato_european_vs_starch_nucleotide_diversity_region.xls | awk 'BEGIN{OFS="\t"}{print $1,$2,$3}' | sort -k1,1 -k2,2n | bedtools merge -d 50000 > 166_potato_european_vs_starch_nucleotide_diversity_region.bed

sort -k1,1 -k2,2n 166_potato_european_vs_starch_nucleotide_diversity_region.bed -o 166_potato_european_vs_starch_nucleotide_diversity_region.bed
bedtools intersect -a 166_potato_european_vs_starch_nucleotide_diversity_region.bed -b DM_v6.1_gene_sort.xls -wb | awk '{print $4,$5,$6}' > 166_potato_european_vs_starch_nucleotide_diversity_region_overlap_gene.xls

rm 166_potato_european_vs_starch_nucleotide_diversity_region_overlap_gene_information.xls
cat 166_potato_european_vs_starch_nucleotide_diversity_region_overlap_gene.xls  | while read l
do
    chr=`echo $l | awk '{print $1}'`
    start=`echo $l | awk '{print $2}'`
    end=`echo $l | awk '{print $3}'`
	grep $chr DM_v6.1_GO_IPR_Swiss_Ara_KEGG.txt | grep $start | grep $end >> 166_potato_european_vs_starch_nucleotide_diversity_region_overlap_gene_information.xls
done
./add_qtl_tpm_2_pai_gene.py reported_qtl_gwas_interval.bed DM_v6.1.tpm.txt 166_potato_european_vs_starch_nucleotide_diversity_region_overlap_gene_information.xls > t && mv t 166_potato_european_vs_starch_nucleotide_diversity_region_overlap_gene_information.xls

### fst per site to MODERATE/HIGH variants, whether in highly diverged regions

vcf=166_potato_Allchr.gatk.filter.SNP.filterMissMAF.vcf

java -jar /media/bulk_01/users/li322/software/snpEff/snpEff.jar -ud 2000 DM_v6.1 $vcf > ${vcf}.snpeff

cat ../head <(egrep 'MODERATE|HIGH' ${vcf}.snpeff) > ${vcf}.MODERATE_HIGH.snpeff

paste <(grep -v '#' ${vcf}.MODERATE_HIGH.snpeff | cut -f '1-5') <(grep -v '#' ${vcf}.MODERATE_HIGH.snpeff|awk '{print $8}'|awk -F 'ANN=' '{print $NF}'|egrep 'MODERATE|HIGH'|awk -F '|' '{print $7}') > ${vcf}.MODERATE_HIGH_gene.xls

rm ${vcf}.MODERATE_HIGH_gene_fst.xls
for chr in {01..12}
do
	grep "chr$chr" ${vcf}.MODERATE_HIGH_gene.xls > ${vcf}.MODERATE_HIGH_gene_chr${chr}.xls
	cat /media/scratchpad_02/li322/work/04_137_atlas_population_analyses/nucleotide_diversity_fst/166_potato_pop_ana_results/chr${chr}/chr${chr}_Fst_each_site.xls3 | ./add_fst.py ${vcf}.MODERATE_HIGH_gene_chr${chr}.xls 166_potato_starch_vs_european_Fst_region.xls chr${chr} >> ${vcf}.MODERATE_HIGH_gene_fst.xls
done

awk '$(NF-1)>0.3' 166_potato_Allchr.gatk.filter.SNP.filterMissMAF.vcf.MODERATE_HIGH_gene_fst.xls | grep -v NA | grep -v Not > 166_potato_Allchr.gatk.filter.SNP.filterMissMAF.vcf.MODERATE_HIGH_gene_fst0.3_inDivergedRegions.xls
#

./calc_allele_freq_gene_func.py pop1_starch.xls pop2_european.xls 166_potato_Allchr.gatk.filter.SNP.filterMissMAF.vcf.MODERATE_HIGH.snpeff 166_potato_Allchr.gatk.filter.SNP.filterMissMAF.vcf.MODERATE_HIGH_gene_fst0.3_inDivergedRegions.xls DM_v6.1_GO_IPR_Swiss_Ara_KEGG.txt reported_qtl_gwas_interval.bed DM_v6.1.tpm.txt > 166_potato_Allchr.gatk.filter.SNP.filterMissMAF.vcf.MODERATE_HIGH_gene_fst0.3_inDivergedRegions_alleleFreq.xls

