#!/usr/bin/env Rscript

library('PopGenome')
argv<-commandArgs(TRUE)

GENOME.class <- readVCF(argv[1],numcols=10000, tid=argv[2], from=1, to=argv[3], approx=FALSE, include.unknown=TRUE, parallel=FALSE, gffpath=argv[4])

#### population assignment, each accession each line
p1 <- as.character(read.table(argv[8])[[1]])
p2 <- as.character(read.table(argv[9])[[1]])

GENOME.class <- set.populations(GENOME.class,list(p1, p2), diploid=FALSE)
slide <- sliding.window.transform(GENOME.class,width=as.numeric(argv[5]),jump=as.numeric(argv[5])/2, type=2, start.pos=1, end.pos=as.numeric(argv[3])+as.numeric(argv[5]))
slide <- diversity.stats(slide)
nucdiv <- slide@nuc.diversity.within
nucdiv <- nucdiv/(100000)

slide <- F_ST.stats(slide, mode="nucleotide")
pairwise.FST <- t(slide@nuc.F_ST.pairwise)

#slide <- neutrality.stats(slide, FAST=FALSE)
#tajimaD <- slide@Tajima.D

write.table(nucdiv,argv[6],row.names = T,col.names = T,quote = FALSE,sep='\t')
write.table(pairwise.FST,argv[7],row.names = T,col.names = T,quote = FALSE,sep='\t')
#write.table(tajimaD,argv[8],row.names = T,col.names = T,quote = FALSE,sep='\t')

