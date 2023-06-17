#!/usr/bin/env Rscript

argv<-commandArgs(TRUE)
#threshold=as.numeric(argv[5]) ### 0.2153, genomewide 5% threshold FST

pdf(argv[2], 10, 3)
par(mar=c(2,4,1,1),oma=c(0,0,0,0),mgp=c(2.5,0.5,0),cex.axis=1.1,las=1,cex.lab=1.2,lend=1)
a <- read.table(argv[1],header=F)

xleft <- min(a[,2])
xright <- max(a[,2])

plot(a$V2,a$V3*1000,type='l',xaxs='i',yaxs='i',cex=0.1,col='#e6c02e',axes=F,main='',xlab='Chromosome',ylab=paste('pai', '× 10-3', sep=' '),xlim=c(xleft, xright),ylim=c(0,15), lwd=2)
par(new=T)

plot(a$V2,a$V4*1000,type='l',xaxs='i',yaxs='i',cex=0.1,col='#c994ec',axes=F,main='',xlab='',ylab='',xlim=c(xleft, xright),ylim=c(0,15), lwd=2)

axis(2,tcl=-0.3,lwd=2,lwd.ticks=2,at=seq(0,15,3),labels=seq(0,15,3))
axis(1,at=seq(xleft, xright, (xright-xleft)/5),labels=seq(xleft/1000000, xright/1000000, (xright-xleft)/5000000),lwd=2,lwd.ticks=2,tcl=-0.3) #### chromosomes

legend('topright',legend=c('Starch', 'Others'),col=c('#e6c02e', '#c994ec'), bty='n',lty=1,lwd=2,cex=1)


