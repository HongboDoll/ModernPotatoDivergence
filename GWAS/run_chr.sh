#!/bin/bash
    . ~/.bashrc
    cd /media/scratchpad_03/li322/work/07_gwas
    conda init
    conda activate R
    

    ./gwaspoly.R 166_potato_Allchr.gatk.filter.SNP.filterMissMAF.genotype_chr12 166_potato_50_phenotype_population_structure.xls 166_potato_50_phenotype_GWASpoly_results.xls_chr12_ 2 4 166_potato_50_phenotype_GWASpoly_qtl_chr12.xls
    

    ./gwaspoly.R 166_potato_Allchr.gatk.filter.SNP.filterMissMAF.genotype_chr09 166_potato_50_phenotype_population_structure.xls 166_potato_50_phenotype_GWASpoly_results.xls_chr09_ 2 4 166_potato_50_phenotype_GWASpoly_qtl_chr09.xls
    

    ./gwaspoly.R 166_potato_Allchr.gatk.filter.SNP.filterMissMAF.genotype_chr07 166_potato_50_phenotype_population_structure.xls 166_potato_50_phenotype_GWASpoly_results.xls_chr07_ 2 4 166_potato_50_phenotype_GWASpoly_qtl_chr07.xls
    

    ./gwaspoly.R 166_potato_Allchr.gatk.filter.SNP.filterMissMAF.genotype_chr05 166_potato_50_phenotype_population_structure.xls 166_potato_50_phenotype_GWASpoly_results.xls_chr05_ 2 4 166_potato_50_phenotype_GWASpoly_qtl_chr05.xls
    

    ./gwaspoly.R 166_potato_Allchr.gatk.filter.SNP.filterMissMAF.genotype_chr04 166_potato_50_phenotype_population_structure.xls 166_potato_50_phenotype_GWASpoly_results.xls_chr04_ 2 4 166_potato_50_phenotype_GWASpoly_qtl_chr04.xls
    

    ./gwaspoly.R 166_potato_Allchr.gatk.filter.SNP.filterMissMAF.genotype_chr01 166_potato_50_phenotype_population_structure.xls 166_potato_50_phenotype_GWASpoly_results.xls_chr01_ 2 4 166_potato_50_phenotype_GWASpoly_qtl_chr01.xls
    
