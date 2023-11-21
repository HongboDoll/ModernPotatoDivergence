#!/bin/bash

ref=DM_v6.1.fa
threads=1

cd /media/scratchpad_02/li322/work/02_137_atlas_potato_reseq

ls mapping/*bam | awk -F '/' '{print $NF}' | awk -F '.' '{print $1}' | while read i
do
    echo """#!/bin/bash
    cd /media/bulk_01/users/li322/work/02_137_atlas_potato_reseq
    gatk --java-options \"-Xmx10G -XX:ParallelGCThreads=$threads\" HaplotypeCaller -R $ref --sample-ploidy 4 --emit-ref-confidence GVCF -I mapping/${i}.sort.markdup.bam -O gvcf/${i}.sample.g.vcf.gz --native-pair-hmm-threads $threads # --sample-ploidy assumes samples are polyploidy
    """ > ${i}_gatk.sh
done

rm sample.gvcf.xls
ls $PWD/gvcf/*.g.vcf.gz | while read i
do
    name=`echo $i | awk -F '/' '{print $NF}' | awk -F '.' '{print $1}'`
    echo -e "$name\t$i" >> sample.gvcf.xls
done

rm vcf.list
grep '>' $ref | grep chr | sed 's/>//g' | while read i
do
    echo """#!/bin/bash
    gatk --java-options \"-Xmx10G -XX:ParallelGCThreads=$threads\" GenomicsDBImport -R $ref -L $i --sample-name-map sample.gvcf.xls --genomicsdb-workspace-path ${i}.gatk.db

    gatk --java-options \"-Xmx10G -XX:ParallelGCThreads=$threads\" GenotypeGVCFs -R $ref --allow-old-rms-mapping-quality-annotation-data --sample-ploidy 2 -V gendb://${i}.gatk.db -O ${i}.gatk.vcf # --sample-ploidy
    echo -e \"$PWD/${i}.gatk.vcf\" >> vcf.list
    """ > ${i}_run.sh && chmod 755 ${i}_run.sh
done


for i in {01..12}
do
echo """#!/bin/bash
	
gatk SelectVariants -select-type SNP --restrict-alleles-to BIALLELIC -V chr${i}.gatk.vcf -O chr${i}.gatk.snp.vcf.gz

gatk --java-options \"-Xmx100G -XX:ParallelGCThreads=13\" VariantFiltration -V ${i}.gatk.snp.vcf.gz \
--filter-expression \"QD < 2.0\" --filter-name \"LowQD\" \
--filter-expression \"MQ < 40.0\" --filter-name \"MQ40.0\" \
--filter-expression \"FS > 60.0\" --filter-name \"FS60.0\" \
--filter-expression \"SOR > 3.0\" --filter-name \"SOR3.0\" \
--filter-expression \"MQRankSum < -12.5\" --filter-name \"MQRankSum-12.5\" \
--filter-expression \"ReadPosRankSum < -8.0\" --filter-name \"ReadPosRankSum-8.0\" \
-G-filter \"GQ < 20.0\" -G-filter-name \"GQ20.0\" \
-O ${i}.gatk.snp.filter.vcf.gz

gatk SelectVariants --max-indel-size 30 --restrict-alleles-to BIALLELIC -select-type INDEL -V ${i}.gatk.vcf -O ${i}.gatk.indel.vcf.gz

gatk --java-options \"-Xmx100G -XX:ParallelGCThreads=13\" VariantFiltration -V ${i}.gatk.indel.vcf.gz \
--filter-expression \"QD < 2.0\" --filter-name \"LowQD\" \
--filter-expression \"MQ < 40.0\" --filter-name \"MQ40.0\" \
--filter-expression \"FS > 200.0\" --filter-name \"FS200\" \
--filter-expression \"SOR > 10.0\" --filter-name \"SOR10\" \
--filter-expression \"MQRankSum < -12.5\" --filter-name \"MQRankSum-12.5\" \
--filter-expression \"ReadPosRankSum < -8.0\" --filter-name \"ReadPosRankSum-8.0\" \
-G-filter \"GQ < 20.0\" -G-filter-name \"GQ20.0\" \
-O ${i}.gatk.indel.filter.vcf.gz

gatk MergeVcfs \
-I ${i}.gatk.snp.filter.vcf.gz \
-I ${i}.gatk.indel.filter.vcf.gz \
-O ${i}.gatk.filter.vcf && rm ${i}.gatk.vcf ${i}.gatk.snp.filter.vcf.gz ${i}.gatk.indel.filter.vcf.gz

cat <(grep '#' ${i}.gatk.filter.vcf | grep -v scaffold) <(grep -v '#' ${i}.gatk.filter.vcf | grep -v scaffold | awk '\$7==\"PASS\"&&\$5!~/,/') > t && mv t ${i}.gatk.filter.vcf

echo -e \"$PWD/${i}.gatk.filter.vcf\" >> vcf.filter.list
""" > ${i}_run.sh && chmod 755 ${i}_run.sh

done

gatk MergeVcfs -I vcf.filter.list -O 166_potato_Allchr.gatk.filter.vcf.gz && rm chr*.gatk.filter.vcf chr*.gatk.vcf

tabix -f -p vcf chr${i}.gatk.filter.vcf.gz

snp=`grep -v '#' 166_potato_Allchr.gatk.filter.vcf | awk 'length($4)==length($5)' | wc -l`
indel=`grep -v '#' 166_potato_Allchr.gatk.filter.vcf | awk 'length($4)!=length($5)' | wc -l`
echo -e "snpNum\tindelNum"
echo -e "${snp}\t${indel}"


#### SNP indel vcf
cat <(grep '#' 166_potato_Allchr.gatk.filter.vcf) <(grep -v '#' 166_potato_Allchr.gatk.filter.vcf | awk 'length($4)==length($5)') > 166_potato_Allchr.gatk.filter.SNP.vcf
cat <(grep '#' 166_potato_Allchr.gatk.filter.vcf) <(grep -v '#' 166_potato_Allchr.gatk.filter.vcf | awk 'length($4)!=length($5)') > 166_potato_Allchr.gatk.filter.indel.vcf

### functional impact
java -jar /media/bulk_01/users/li322/software/snpEff/snpEff.jar -ud 2000 DM_v6.1 166_potato_Allchr.gatk.filter.vcf > 166_potato_Allchr.gatk.filter.vcf.snpeff
gzip 166_potato_Allchr.gatk.filter.vcf.snpeff
java -jar /media/bulk_01/users/li322/software/snpEff/snpEff.jar -ud 2000 DM_v6.1 166_potato_Allchr.gatk.filter.SNP.vcf > 166_potato_Allchr.gatk.filter.SNP.vcf.snpeff && mv snpEff_summary.html snp_vcf_snpEff_summary.html
gzip 166_potato_Allchr.gatk.filter.SNP.vcf.snpeff
java -jar /media/bulk_01/users/li322/software/snpEff/snpEff.jar -ud 2000 DM_v6.1 166_potato_Allchr.gatk.filter.indel.vcf > 166_potato_Allchr.gatk.filter.indel.vcf.snpeff && mv snpEff_summary.html indel_vcf_snpEff_summary.html
gzip 166_potato_Allchr.gatk.filter.indel.vcf.snpeff


### observed heterozygosity: at a given site, number of homozygotes / number of individuals without missing genotype
./output_hetero_variants.py 166_potato_Allchr.gatk.filter.vcf 0.5 > 166_potato_Allchr.gatk.filter.heterozygous0.5.vcf
grep -v '#' 166_potato_Allchr.gatk.filter.SNP.heterozygous0.5.vcf | wc -l > 166_potato_Allchr.gatk.filter.SNP.heterozygous0.5.vcf.line

zcat 166_potato_Allchr.gatk.filter.vcf.gz | ./calc_observed_expected_hetero.py > 166_potato_Allchr.gatk.filter.vcf.observed_expected_hetero.xls

mean_ob_hetero=`awk '{i+=$(NF-1)}END{print i/NR}' 166_potato_Allchr.gatk.filter.vcf.observed_expected_hetero.xls`
mean_exp_hetero=`awk '{i+=$NF}END{print i/NR}' 166_potato_Allchr.gatk.filter.vcf.observed_expected_hetero.xls`

echo $mean_ob_hetero $mean_exp_hetero

### SNP indel density
grep -v '#' 166_potato_Allchr.gatk.filter.vcf | awk 'BEGIN{OFS="\t"}{print $1"_"$2,$1,$2}' > 166_potato_Allchr.gatk.filter_4_CMplot.xls

./plot_snp_density.R 166_potato_Allchr.gatk.filter_4_CMplot.xls 166_potato_Allchr_SNP_density500k

