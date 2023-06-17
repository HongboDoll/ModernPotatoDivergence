#!/bin/bash

cd /media/scratchpad_02/li322/work/04_137_atlas_population_analyses/nucleotide_diversity_fst/fst_permutation

window=100000   # argv[5]
slide=50000   # argv[5]/2, half of the window size
gff=DM_v6.1.gff3   # argv[4]

#for n in `seq 1 10`
#do
#shuf -n23 pop.xls > pop1_${n}.xls
#shuf -n102 pop.xls > pop2_${n}.xls
#
#rm -rf 166_potato_pop_ana_results${n}; mkdir 166_potato_pop_ana_results${n}
#grep '>chr' DM_v6.1.fa | sed 's/>//g' | while read i
#do
#       mkdir -p 166_potato_pop_ana_results${n}/${i}
#       end=`out_len.py DM_v6.1.fa | grep $i | awk '{print $2}'`  # argv[3]
#       ref=$i   # argv[2]
#
#       nuc_div=${i}_nucleotide_diversity.xls   # argv[6]
#       fst=${i}_Fst.xls   # argv[7]
#       vcf=$PWD/filter_vcf/${i}/166_potato_Allchr.gatk.filter.SNP.filterMissMAF.${i}.vcf.gz   # argv[1]
#	   echo -e """#!/bin/bash
#
#	   cd /media/scratchpad_02/li322/work/04_137_atlas_population_analyses/nucleotide_diversity_fst/fst_permutation/166_potato_pop_ana_results${n}/${i}/
#
#       ../../popgenome.R $vcf $ref $end ../../$gff $window $nuc_div $fst ../../pop1_${n}.xls ../../pop2_${n}.xls
#	   """ > ${i}_pop_ana.sh && chmod 755 ${i}_pop_ana.sh
#		qsub -cwd -N ${i}_pai_fst -pe sharedmem 1 -q seq.q -j y -V -b y -S /bin/bash "$PWD/${i}_pop_ana.sh"
#		
#done
#	sleep 30m
#done


window=100000
slide=50000

for n in `seq 1 10`
do

grep '>chr' DM_v6.1.fa | sed 's/>//g' | while read i
do

       ref=$i   # argv[2]
       nuc_div=166_potato_pop_ana_results${n}/${i}/${i}_nucleotide_diversity.xls   # argv[6]
       fst=166_potato_pop_ana_results${n}/${i}/${i}_Fst.xls   # argv[7]

#       awk 'NR!=1' 166_potato_pop_ana_results${n}/${i}/${i}_nucleotide_diversity.xls | awk -F '\t' 'BEGIN{OFS="\t"}{print "'"$i"'",($1-1)*"'"$slide"'",($1-1)*"'"$slide"'"+"'"$window"'",$2,$3}' > t && mv t 166_potato_pop_ana_results${n}/${i}/${i}_nucleotide_diversity.xls
#       awk 'NR!=1' 166_potato_pop_ana_results${n}/${i}/${i}_Fst.xls | awk -F '\t' 'BEGIN{OFS="\t"}{print "'"$i"'",($1-1)*"'"$slide"'",($1-1)*"'"$slide"'"+"'"$window"'",$2}' > t && mv t 166_potato_pop_ana_results${n}/${i}/${i}_Fst.xls

done

#cat <(echo -e "Chr\tPos\tPos2\tPop1\tPop2") <(cat 166_potato_pop_ana_results${n}/*/*_nucleotide_diversity.xls) > 166_potato_pop1_pop2_nucleotide_diversity_${n}.xls
#cat <(echo -e "Chr\tPos\tPos2\tPop1_vs_Pop2") <(cat 166_potato_pop_ana_results${n}/*/*_Fst.xls) > 166_potato_pop1_pop2_Fst_${n}.xls

cat <(cut -f '1,2,4' 166_potato_pop1_pop2_Fst_${n}.xls | head -1) <(cut -f '1,2,4' 166_potato_pop1_pop2_Fst_${n}.xls | awk 'NR!=1' | awk 'BEGIN{OFS="\t"}{if($3<0){print $1,$2,0}else{print $1,$2,$3}}') > t
./plot_fst.R t 166_potato_pop1_pop2_Fst_${n}.pdf 12 chr 0.2135

cat <(echo -e "Chr\tPos\tPai") <(awk 'NR!=1' 166_potato_pop1_pop2_nucleotide_diversity_${n}.xls |awk '$5>0.0002&&$4>0.0001'|awk '{if($5!=0){print $1,$2,$4/$5}else{print $1,$2,0}}') > t
./plot_pai.R t 166_potato_pop1_pop2_nucleotide_diversity_${n}.pdf 12 chr 1.84654

awk 'NR!=1' 166_potato_pop1_pop2_nucleotide_diversity_${n}.xls |awk '{if($5>0.0002&&$4>0.0001){print $0}else{print $1,$2,$3,0,0}}'|awk '{if($4!=0){print $1,$2,$5/$4}else{print $1,$2,0}}' | awk '{if($3>1.84654){print 1}else{print 0}}' > pai_permutation_${n}.xls

done

awk 'NR!=1' 166_potato_3pop_nucleotide_diversity.xls |awk '{if($5>0.0002&&$4>0.0001){print $0}else{print $1,$2,$3,0,0}}'|awk '{if($4!=0){print $1,$2,$5/$4}else{print $1,$2,0}}' | awk '{if($3>1.84654){print 1}else{print 0}}' > european_starch_pai.xls

cat head <(paste european_starch_pai.xls pai_permutation_*.xls) > european_starch_permutation_1-10_pai.xls



