#!/usr/bin/env Rscript

argv<-commandArgs(TRUE)
d <- read.table(argv[1],header=T)

e <- p.adjust(d$fisher_p_value,method='fdr',length(d$fisher_p_value))

write.table(e,argv[2],row.names = FALSE,col.names = FALSE,quote = FALSE,sep='\n')

	

