# 15:5
d <- read.table('whole_genome_combined_D_Z_order.xls', header=F, sep='\t')
row.names(d) <- d$V1
par(mgp=c(3,0.8,0),cex.lab=1.5,cex.axis=1.3,las=1,oma=c(0,0,0,0), mar=c(9.5,6.2,1,1), lend=1)

### 10 colors
colors <- rev(colorRampPalette(c("#88419d", "#e0ecf4"))(10))
plot(d$V2, pch=21, col = colors[rank(-(d$V3))], bg = colors[rank(-(d$V3))],cex = 0, axes=F, xaxs='i',yaxs='i',ylim=c(0,0.3), xlim=c(0, length(d[,1])+1),ylab='Patterson\'s D', xlab='')

### divide d$V3 (z-score) to 10 intervals
interval <- (ceiling(max(d[,3]))-floor(min(d[,3])))/10
range <- c()
for (i in 1:10){
    range <- c(range, floor(min(d[,3]))+interval*(i-1))
    range <- c(range, floor(min(d[,3]))+interval*(i))
}

for (i in 1:length(d[,1])){
    f <- 0
    for (n in 1:10){
        if (d[i,3] >= range[2*n-1] && d[i,3] < range[2*n]){
            f = f + n
        }
    }
    points(i, d[i,2], cex=1.6, pch=21,  col=F, bg=colors[f])
}

axis(1,at=seq(0, length(d[,1]), 1),labels=NA,tcl=-0.5,lwd=2,lwd.ticks=2)
text(x=seq(1, length(d[,1]), 1),y=-0.01, srt = 45, adj = 1, labels = row.names(d),xpd = TRUE,cex=1.3)

axis(2,at=seq(0, 0.3, 0.1),labels=seq(0, 0.3, 0.1),tcl=-0.5,lwd=2,lwd.ticks=2)

range