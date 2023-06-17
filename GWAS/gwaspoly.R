#!/usr/bin/env Rscript

argv<-commandArgs(TRUE)

library(GWASpoly)

pheno <- read.table(argv[2], header=T)
data <- read.GWASpoly(ploidy=4,geno.file=argv[1],pheno.file=argv[2],format="ACGT",n.traits=(length(colnames(pheno))-1-as.numeric(argv[5])),delim="\t")

data_loco <- set.K(data, n.core=as.numeric(argv[4]), LOCO=TRUE)
params <- set.params(fixed=c("Grp1","Grp2","Grp3","Grp4"),n.PC=5,fixed.type=rep("numeric",4), MAF=0.05)

data_loco.scan <- GWASpoly(data=data_loco,models=c("additive","1-dom"),params=params,n.core=as.numeric(argv[4]))

data_loco.scan.threshold <- set.threshold(data_loco.scan,method="Bonferroni",level=0.05)
qtl <- get.QTL(data_loco.scan.threshold,bp.window=10000)

write.table(qtl, argv[6],quote = FALSE,sep='\t')

trait_list <- as.vector(colnames(pheno))[2:(length(colnames(pheno))-as.numeric(argv[5]))]
for (i in 1:length(trait_list)){
	write.GWASpoly(data_loco.scan.threshold, trait=trait_list[i],filename=paste(argv[3],trait_list[i], sep=''), what = "scores", delim = "\t")	
}

