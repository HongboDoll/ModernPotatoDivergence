#!/usr/bin/env Rscript

argv<-commandArgs(TRUE)

pdf(argv[4], 5,5)
par(mar=c(4,4,1,1),oma=c(0,0,0,0),mgp=c(2.5,0.8,0),cex.axis=1.1,las=1,cex.lab=1.2,lend=1)

read.table(argv[1])->Epop1_starch;
plot(Epop1_starch[,1]/1000,Epop1_starch[,2],type="l",col="#e6c02e",main="",xlab="Pairwise distance (kb)",xlim=c(0,50),ylim=c(0.1,0.5),ylab=expression(r^{2}),bty="n",lwd=2, axes=F)
read.table(argv[2])->Epop2_european;
lines(Epop2_european[,1]/1000,Epop2_european[,2],col="#46dff0",lwd=2)
read.table(argv[3])->Epop3_american;
lines(Epop3_american[,1]/1000,Epop3_american[,2],col="#e90c59",lwd=2)

axis(2,tcl=-0.5,lwd=2,lwd.ticks=2,at=seq(0.1,0.5,0.1),labels=seq(0.1,0.5,0.1))
axis(1,at=seq(0, 50, 10),labels=seq(0, 50, 10),lwd=2,lwd.ticks=2,tcl=-0.5) #### chromosomes

#legend("topright",c("pop1_starch","pop2_european","pop3_american"),col=c("red","black","blue"),cex=1,lty=c(1,1,1),bty="n",lwd=2);
legend('topright',legend=c('Starch', 'European', 'American'),col=c('#e6c02e', '#46dff0', '#e90c59'), bty='n',lty=1,lwd=2,cex=1.2)

dev.off()
#png("166_potato_filter_population_LDdecay.png")
#read.table("166_potato_filter_population_LDdecay.pop1_starch")->Epop1_starch;
#plot(Epop1_starch[,1]/1000,Epop1_starch[,2],type="l",col="red",main="LD decay",xlab="Distance(Kb)",xlim=c(0,500),ylim=c(0,0.459050128956662),ylab=expression(r^{2}),bty="n",lwd=2)
#read.table("166_potato_filter_population_LDdecay.pop2_european")->Epop2_european;
#lines(Epop2_european[,1]/1000,Epop2_european[,2],col="black",lwd=2)
#read.table("166_potato_filter_population_LDdecay.pop3_american")->Epop3_american;
#lines(Epop3_american[,1]/1000,Epop3_american[,2],col="blue",lwd=2)
#
#legend("topright",c("pop1_starch","pop2_european","pop3_american"),col=c("red","black","blue"),cex=1,lty=c(1,1,1),bty="n",lwd=2);
#dev.off()
#
