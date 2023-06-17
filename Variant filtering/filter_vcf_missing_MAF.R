#!/usr/bin/env Rscript

#library(ape)
library(PopGenome)
argv<-commandArgs(TRUE)

# argv[1]: chromosome folder vcfs/chr01 chr02 .. chr12
# argv[2]: missing rate threshold 0.2
# argv[3]: MAF threshold 0.05

GENOME.class <- readData(argv[1], format = 'VCF', include.unknown=TRUE)

GENOME.class <- set.filter(GENOME.class, minor.freqs=TRUE,maf.lower.bound=as.numeric(argv[3]),maf.upper.bound=1)

GENOME.class <- set.filter(GENOME.class, missing.freqs=TRUE,miss.lower.bound=0, miss.upper.bound=as.numeric(argv[2]))

l <- c()
for (n in 1:length(GENOME.class@region.data@included[[1]])) {
    if (as.character(GENOME.class@region.data@included[[1]][n]) == "TRUE"){
        l <- c(l, paste(argv[1], GENOME.class@region.data@biallelic.sites[[1]][n], sep='\t'))
    }
}

write.table(l,argv[4],row.names = FALSE,col.names = FALSE,quote = FALSE,sep='')



