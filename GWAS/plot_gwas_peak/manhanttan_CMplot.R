#!/usr/bin/env Rscript

argv<-commandArgs(TRUE)

library('CMplot')

a <- read.table(argv[1], header=F)

CMplot(a,plot.type = c('m','q'),file='pdf',file.name=argv[2],col=c(rgb(129,196,240,max=255),rgb(235,70,144,max=255)),threshold = as.numeric(argv[3]),threshold.col='black',threshold.lty = 2,threshold.lwd = 1, amplify = F,signal.cex = c(1,1), signal.pch = c(20,20),signal.col = NULL, dpi=326, cex.axis=1.3,cex.lab=1.5,cex=0.8,conf.int=T,conf.int.col='grey',box=F,main='',ylab=expression(-log[10](italic(p))))

