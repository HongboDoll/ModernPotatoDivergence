#!/usr/bin/env Rscript

argv<-commandArgs(TRUE)

pdf(argv[2], 5, 5)
dd <- read.table(argv[1],header=F)
par(mgp=c(2.4,0.9,0),las=1,mar=c(3.4,3.6,1,1),oma=c(1,1,1,1),cex.lab=1.2,cex.axis=1.1,lend=1)

sdpc1 <- sd(dd[,2]) * sd(dd[,2])
sdpc2 <- sd(dd[,3]) * sd(dd[,3])
sdsum <- 0
for (i in 2:length(dd[1,])){
    sdpc <- sd(dd[,i]) * sd(dd[,i])
    sdsum <- sdsum + sdpc
}
pc1_variance_explained <- round(sdpc1/sdsum * 100, 2)
pc2_variance_explained <- round(sdpc2/sdsum * 100, 2)


plot(dd$V2,dd$V3,pch=rep(20,166),col=c(rep("#2cbefb",74),rep("black",12), rep("#c994ec",47), rep("#e6c02e",23),rep("#83d28f",10)),xlab=paste("PC1 (", pc1_variance_explained, "% of variance)", sep=''),ylab=paste("PC2 (", pc2_variance_explained, "% of variance)", sep=''),bty="n",tcl=-0.2,xaxs='i',yaxs='i',axes=F,xlim=c(-40,60),ylim=c(-40,40),cex=1.2) # fresh, others, Processing, starch, outgroup
axis(2,at=seq(-40,40,20),labels=seq(-40,40,20),col.axis='black',col.lab='black',lwd.ticks=2,lwd=2)
axis(1,at=seq(-40,60,20),labels=seq(-40,60,20),col.axis='black',col.lab='black',lwd.ticks=2,lwd=2)

legend(25,43,c("Fresh","Processing","Starch","Others","Outgroup"),pch=rep(19,5),col=c("#2cbefb","#c994ec","#e6c02e","black","#83d28f"),bty='n',cex=1.1)

