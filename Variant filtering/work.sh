#!/bin/bash

#cd /media/scratchpad_02/li322/work/02_137_atlas_potato_reseq
#bgzip -d 166_potato_Allchr.gatk.filter.SNP.vcf.gz
#cd -
#ln -s /media/scratchpad_02/li322/work/02_137_atlas_potato_reseq/166_potato_Allchr.gatk.filter.SNP.vcf


cd /media/scratchpad_02/li322/work/04_137_atlas_population_analyses

iTools Fatools  getCdsPep -Ref DM_v6.1.fa -Gff DM_v6.1.gff3 -4DSite -OutPut DM_v6.1_4dSNP
gzip -d DM_v6.1_4dSNP.4Dsite.gz

rm -rf vcfs/; mkdir vcfs
for i in {01..12}
do
    mkdir -p vcfs/chr${i}
	echo -e """#!/bin/bash
	cd /media/scratchpad_02/li322/work/04_137_atlas_population_analyses
    cat <(grep '#' 166_potato_Allchr.gatk.filter.SNP.vcf) <(grep -v '#' 166_potato_Allchr.gatk.filter.SNP.vcf | grep chr${i}) > vcfs/chr${i}/166_potato_Allchr.gatk.filter.SNP.chr${i}.vcf
	""" > vcf_chr${i}.sh && chmod 755 vcf_chr${i}.sh
	qsub -cwd -N chr${i} -q stat.q -j y -V -b y -S /bin/bash "//media/scratchpad_02/li322/work/04_137_atlas_population_analyses/vcf_chr${i}.sh"
done

head -500 vcfs/chr01/166_potato_Allchr.gatk.filter.SNP.chr01.vcf | grep '#' > head_vcf
rm -rf split_vcf
for i in {01..12}
do
    mkdir -p split_vcf/chr${i}
    split -a 3 -d -l 200000 <(grep -v '#' vcfs/chr${i}/166_potato_Allchr.gatk.filter.SNP.chr${i}.vcf) split_vcf/chr${i}/166_potato_Allchr.gatk.filter.SNP.chr${i}.vcf_
done

rm -rf filter_vcf 4d_vcf; mkdir filter_vcf 4d_vcf
rm -rf filter_retained_sites; mkdir filter_retained_sites
for i in {01..12}
do
    ls split_vcf/chr${i} | while read v
    do
        num=`echo $v | awk -F '_' '{print $NF}'`
		cat head_vcf split_vcf/chr${i}/$v > tmp && mv tmp split_vcf/chr${i}/$v
		mkdir -p split_vcf/chr${i}/${num} && mv split_vcf/chr${i}/$v split_vcf/chr${i}/${num}
        ./filter_vcf_missing_MAF.R split_vcf/chr${i}/${num} 0.2 0.05 filter_retained_sites/filter_retained_sites_chr${i}_${num}
    done
    cat filter_retained_sites/filter_retained_sites_chr${i}_* | awk 'BEGIN{OFS="\t"}{print "chr""'"$i"'",$2}' > filter_retained_sites_chr${i}
	mkdir -p filter_vcf/chr${i}
#	mkdir -p 4d_vcf/chr${i}
    ./retain_vcf_sites.py vcfs/chr${i}/166_potato_Allchr.gatk.filter.SNP.chr${i}.vcf filter_retained_sites_chr${i} > filter_vcf/chr${i}/166_potato_Allchr.gatk.filter.SNP.filterMissMAF.chr${i}.vcf
#    ./retain_vcf_sites.py filter_vcf/chr${i}/166_potato_Allchr.gatk.filter.SNP.filterMissMAF.chr${i}.vcf DM_v6.1_4dSNP.4Dsite > 4d_vcf/chr${i}/166_potato_Allchr.gatk.filter.SNP.filterMissMAF.4d.chr${i}.vcf
done


cat head_vcf <(cat filter_vcf/chr*/166_potato_Allchr.gatk.filter.SNP.filterMissMAF.chr*.vcf | grep -v '#' | sort -k1,1 -k2,2n) > 166_potato_Allchr.gatk.filter.SNP.filterMissMAF.vcf
#
##cat head_vcf <(cat 4d_vcf/chr*/166_potato_Allchr.gatk.filter.SNP.filterMissMAF.4d.chr*.vcf  | grep -v '#' | sort -k1,1 -k2,2n) > 166_potato_Allchr.gatk.filter.SNP.filterMissMAF.4d.vcf
#
