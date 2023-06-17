#!/bin/bash

cd /media/scratchpad_02/li322/work/04_137_atlas_population_analyses/phylogeny

vcf=166_potato_Allchr.gatk.filter.SNP.filterMissMAF.4d.vcf

#./convert_4d_vcf_2_allele_freq_4_phylip.py ${vcf} corresondance > ${vcf}.4_phylip

### calculate genetic distance among accessions based on allele frequency (dependent of ploidy levels) and construct phylogenetic relationships using phylip (NJ method)
rm outfile
#seqboot < seqboot.params
#mv outfile bootstrap1000.out
#gendist < gendist.params
#mv outfile gendist1000.out
neighbor < neighbor.params
mv outtree intree
consense < consense.params
mv outtree consensus_tree.tre
