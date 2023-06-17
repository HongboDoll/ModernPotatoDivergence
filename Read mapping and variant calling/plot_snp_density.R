#!/usr/bin/env Rscript
argv<-commandArgs(TRUE) 

library('CMplot')
a <- read.table(argv[1])

CMplot(a, plot.type="d",  bin.size=5e5, chr.den.col=c(rgb(129,196,240,max=255),rgb(255,47,146,max=255)),
       file="jpg", memo=argv[2], dpi=326, file.output=TRUE, verbose=TRUE, main='')

