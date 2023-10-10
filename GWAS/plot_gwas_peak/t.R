#!/usr/bin/env Rscript

argv<-commandArgs(TRUE)
a <- read.table(argv[1],fill=T)

write.table(t(a),file=argv[2],row.names = F,col.names = FALSE,quote = FALSE,sep='\t')



