#!/usr/bin/env Rscript

argv<-commandArgs(TRUE)

### import data
library(adegenet)
a <- read.table(argv[1],header=T) # 
obj <- df2genind(a, ploidy=4, sep="/",NA.char='NA',pop=c(rep('Indian', 30),rep('EA',84))) # pop means assigned sub-populations according to prior knowledge

### PCA
pca <- tab(obj, freq = TRUE, NA.method = "mean")
pca.out <- dudi.pca(pca, scale = FALSE, scannf = FALSE, nf = 20)
write.table(pca.out$li,argv[2],row.names = T,col.names = F,quote = FALSE,sep=' ') # output H-W test p-value 

