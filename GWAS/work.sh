#!/bin/bash


cd /media/scratchpad_03/li322/work/07_gwas

######  GEC threshold

cd /media/scratchpad_02/li322/work/07_gwas

./convert_vcf_2_gwaspoly_genotype.py 166_potato_Allchr.gatk.filter.SNP.filterMissMAF.vcf > 166_potato_Allchr.gatk.filter.SNP.filterMissMAF.genotype

cp 166_potato_Allchr.gatk.filter.SNP.filterMissMAF.genotype /media/bulk_01/users/li322/work/07_gwas

plink --vcf 166_potato_Allchr.gatk.filter.SNP.filterMissMAF.diploid_genotype.vcf --recode12 --allow-extra-chr --geno 0.2 --maf 0.05 --allow-no-sex --biallelic-only --out 166_potato

plink --file 166_potato --allow-extra-chr --make-bed --out 166_potato

java -jar /media/bulk_01/users/li322/software/gec/gec.jar -Xmxlg --effect-number --plink-binary ./166_potato --genome --no-web --out 166_potato.gec

threshold=`grep -v "Observed_Number" 166_potato.gec.sum | awk '{print $4}'`


paste 166_potato_50_phenotype.xls <(cat <(echo -e "name\tGrp1\tGrp2\tGrp3\tGrp4") 166_potato_structure_K4.xls | cut -f '2-99') > t && mv t 166_potato_50_phenotype_population_structure.xls

./convert_vcf_2_gwaspoly_genotype.py 166_potato_Allchr.gatk.filter.SNP.filterMissMAF.vcf > 166_potato_Allchr.gatk.filter.SNP.filterMissMAF.genotype

cp 166_potato_Allchr.gatk.filter.SNP.filterMissMAF.genotype /media/bulk_01/users/li322/work/07_gwas

###### only can be run on seq.long.q
n_populations=4
n_cores=2
#for i in {01..12}
rm run_chr.sh
#   cat <(head -1 166_potato_Allchr.gatk.filter.SNP.filterMissMAF.genotype) <(grep chr${i} 166_potato_Allchr.gatk.filter.SNP.filterMissMAF.genotype) > 166_potato_Allchr.gatk.filter.SNP.filterMissMAF.genotype_chr${i}

    #nohup  ./gwaspoly.R 166_potato_Allchr.gatk.filter.SNP.filterMissMAF.genotype_chr${i} 166_potato_50_phenotype_population_structure.xls 166_potato_50_phenotype_GWASpoly_results.xls_chr${i}_ $n_cores $n_populations 166_potato_50_phenotype_GWASpoly_qtl_chr${i}.xls &
    echo -e """#!/bin/bash
    . ~/.bashrc
    cd /media/scratchpad_03/li322/work/07_gwas
    conda init
    conda activate R
    """ > run_chr.sh
for i in 12 09 07 05 04 01
do
    echo -e """
    ./gwaspoly.R 166_potato_Allchr.gatk.filter.SNP.filterMissMAF.genotype_chr${i} 166_potato_50_phenotype_population_structure.xls 166_potato_50_phenotype_GWASpoly_results.xls_chr${i}_ $n_cores $n_populations 166_potato_50_phenotype_GWASpoly_qtl_chr${i}.xls
    """ >> run_chr.sh && chmod 755 run_chr.sh
    #qsub -cwd -N gwas -q seq.q -j y -V -b y -pe sharedmem 2 -S /bin/bash "$PWD/run_chr${i}.sh"


done

    qsub -cwd -N gwas2 -q seq.long.q -j y -V -b y -pe sharedmem 24 -S /bin/bash "$PWD/run_chr.sh"

####### plot

rm -rf plot plot.sh; mkdir plot
cat 50_phenotype_list.xls | while read trait
do
    cat 166_potato_50_phenotype_GWASpoly_results.xls_chr*_${trait} >  166_potato_50_phenotype_GWASpoly_results.xls_Allchr_${trait}
#   for m in `seq 6 8`
#   do
        awk 'NR!=1' 166_potato_50_phenotype_GWASpoly_results.xls_Allchr_${trait} | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,10^(-$6)}' > plot/166_potato_50_phenotype_GWASpoly_results.xls_Allchr_${trait}_additive_4_plot
        awk 'NR!=1' 166_potato_50_phenotype_GWASpoly_results.xls_Allchr_${trait} | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,10^(-$7)}' > plot/166_potato_50_phenotype_GWASpoly_results.xls_Allchr_${trait}_1-dom-alt_4_plot
        awk 'NR!=1' 166_potato_50_phenotype_GWASpoly_results.xls_Allchr_${trait} | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,10^(-$8)}' > plot/166_potato_50_phenotype_GWASpoly_results.xls_Allchr_${trait}_1-dom-ref_4_plot
echo -e """

        ./manhanttan_CMplot.R plot/166_potato_50_phenotype_GWASpoly_results.xls_Allchr_${trait}_additive_4_plot ${trait}_additive_gwas_plot.jpg $threshold
        ./manhanttan_CMplot.R plot/166_potato_50_phenotype_GWASpoly_results.xls_Allchr_${trait}_1-dom-alt_4_plot ${trait}_1-dom-alt_gwas_plot.jpg $threshold
        ./manhanttan_CMplot.R plot/166_potato_50_phenotype_GWASpoly_results.xls_Allchr_${trait}_1-dom-ref_4_plot ${trait}_1-dom-ref_gwas_plot.jpg $threshold
""" >> plot.sh
#   done
done

split -a 2 -d -l 14 plot.sh plot.sh_

ls plot.sh_* | while read i
do
    cat <(echo -e "#!/bin/bash") <(echo -e "cd /media/scratchpad_02/li322/work/07_gwas") $i > tmp && mv tmp $i && chmod 755 $i
    qsub -cwd -N $i -q seq.q -j y -V -b y -pe sharedmem 1 -S /bin/bash "$PWD/$i"
done

for i in `seq 1 21`
do
    qsub -cwd -N plot.sh_${i} -q seq.q -j y -V -b y -pe sharedmem 1 -S /bin/bash "$PWD/plot.sh_${i}"
done

ls plot/*4_plot| while read i
do
    sort -k4,4g $i > ${i}_sort
done

#mv *jpg plot

#cat 166_potato_50_phenotype_GWASpoly_qtl_chr*.xls | sort -k2,2 -k6,6 -k7,7n > 166_potato_50_phenotype_GWASpoly_qtl_Allchr.xls

###### peak region genes

paste <(echo -e "Trait") <(head -1 DM_v6.1_GO_IPR_Swiss_Ara_KEGG.txt) > peak_region_gene.xls
cat DM_v6.1_GO_IPR_Swiss_Ara_KEGG.txt | grep chr06 | awk '$2>=50131166&&$2<=50725690' | awk '{print "Blackspot bruise\t"$0}' >> peak_region_gene.xls

cat DM_v6.1_GO_IPR_Swiss_Ara_KEGG.txt | grep chr05 | awk '$2>=4395833&&$2<=4895833' | awk '{print "Maturity\t"$0}' >> peak_region_gene.xls
cat DM_v6.1_GO_IPR_Swiss_Ara_KEGG.txt | grep chr09 | awk '$2>=59051717&&$2<=59962076' | awk '{print "Maturity\t"$0}' >> peak_region_gene.xls

cat DM_v6.1_GO_IPR_Swiss_Ara_KEGG.txt | grep chr03 | awk '$2>=42627114&&$2<=43127114' | awk '{print "Tuber flesh color\t"$0}' >> peak_region_gene.xls
cat DM_v6.1_GO_IPR_Swiss_Ara_KEGG.txt | grep chr06 | awk '$2>=3865309&&$2<=4365309' | awk '{print "Tuber flesh color\t"$0}' >> peak_region_gene.xls
cat DM_v6.1_GO_IPR_Swiss_Ara_KEGG.txt | grep chr07 | awk '$2>=55883697&&$2<=56383697' | awk '{print "Tuber flesh color\t"$0}' >> peak_region_gene.xls
cat DM_v6.1_GO_IPR_Swiss_Ara_KEGG.txt | grep chr11 | awk '$2>=43283768&&$2<=43783768' | awk '{print "Tuber flesh color\t"$0}' >> peak_region_gene.xls

cat DM_v6.1_GO_IPR_Swiss_Ara_KEGG.txt | grep chr05 | awk '$2>=52666343&&$2<=53166343' | awk '{print "Tuber skin color\t"$0}' >> peak_region_gene.xls
cat DM_v6.1_GO_IPR_Swiss_Ara_KEGG.txt | grep chr06 | awk '$2>=36286447&&$2<=37074161' | awk '{print "Tuber skin color\t"$0}' >> peak_region_gene.xls
cat DM_v6.1_GO_IPR_Swiss_Ara_KEGG.txt | grep chr10 | awk '$2>=51151931&&$2<=53872620' | awk '{print "Tuber skin color\t"$0}' >> peak_region_gene.xls

cat DM_v6.1_GO_IPR_Swiss_Ara_KEGG.txt | grep chr01 | awk '$2>=59627080&&$2<=60349412' | awk '{print "Resistance to common scab\t"$0}' >> peak_region_gene.xls

###### leading snps

head -1000 166_potato_Allchr.gatk.filter.SNP.filterMissMAF.vcf | grep '#' > vcf_head
rm leading_snp.vcf
cat leading_snp.xls | while read i
do
    grep $i 166_potato_Allchr.gatk.filter.SNP.filterMissMAF.vcf >> leading_snp.vcf
done

cat vcf_head leading_snp.vcf > t && mv t leading_snp.vcf
