######### device size ratio
library('bezier')

a <- read.table('blue_gene_coordinate.xls',header=F)

par(oma=c(0,0,0,0),mar=c(2,4,0,0))

left_margin <- a[1,3] # coordinates of 9930 reference #min(a$V3,b$V3,c$V3,d$V3,e$V3,f$V3), here a[1,3] refer to the left start coordinate of the first gene of 9930

spe_num <- max(as.vector(a[,length(a[1,])-1]))

width <- 0.9
sv_width <- 1.6

#### refine the coordinates of the other species according to Heinz1706 ref
left_l <- c()
for (i in 1:spe_num){
    a1 <- subset(a,V6==i)
    if (a1[1,3] > left_margin){
        left <- min(a1[,3])
        left_l <- c(left_l, left)
        for (n in 1:length(a1[,1])){
            a1[n,3] <- a1[n,3] - (left - left_margin)
            a1[n,4] <- a1[n,4] - (left - left_margin)
        } 
    } else {
        left <- min(a1[,3])
        left_l <- c(left_l, left)
        for (n in 1:length(a1[,1])){
            a1[n,3] <- a1[n,3] + (left_margin - left)
            a1[n,4] <- a1[n,4] + (left_margin - left)	
        }
    }
    a <- rbind(a1, subset(a,V6!=i))
}


#### plot frame and line 
right_margin <- max(a$V4)
inter <- (right_margin - left_margin)/30

#### move all these sgements to middle
xl <- c()
for (i in 1:spe_num){
    a1 <- subset(a,V6==i)
    x <- (right_margin - max(a1[,4]) - min(a1[,3]) + left_margin)/2
    xl <- c(xl, x)
    for (n in 1:length(a1[,1])){
        a1[n,3] <- a1[n,3] + x
        a1[n,4] <- a1[n,4] + x
    }
    a <- rbind(a1, subset(a,V6!=i))
}

#### plot line

plot(0,0,type='n',axes=F,main='',xlab='',ylab='',xlim=c(left_margin - inter, right_margin + inter), ylim=c(1*10-5,spe_num*10+5)) # ylim: min-5=1*10-5, max+5=(spe_num*10)+5

#### all sgenments are in equal length 

#for (i in 1:spe_num){
#    segments(left_margin - inter,i*10,right_margin + inter,i*10,lwd=5,lend=0, #col=rgb(128,128,128,max=255))
#}

#### sgenments lengths are associated with actual gene cluster length

for (i in 1:spe_num){
    a1 <- subset(a,V6==i)
    left <- min(a1[,3])
    right <- max(a1[,4])
    segments(left - inter,(spe_num+1-i)*10,right + inter,(spe_num+1-i)*10,lwd=3.5,lend=1, col=rgb(128,128,128,max=255))
}


#### plot gene body (according to plus and minus strand)
for (n in 1:spe_num){
    a1 <- subset(a,V6==n)
    nn <- n + (spe_num-(2*(n-1)+1))
    for (i in 1:length(a1[,1])){
        if (a1[i,7] == "yellow"){
            gene_col = '#FEDD7F'
        } else if (a1[i,7] == "red") {
            gene_col = 'red'
        } else {
            gene_col = '#eb4690'
        }
        if (a1[i,5] == "+"){
            polygon(c(a1[i,3],a1[i,3],a1[i,3]+0.85*(a1[i,4]-a1[i,3]),a1[i,4],a1[i,3]+0.85*(a1[i,4]-a1[i,3]),a1[i,3]), c(nn*10+width,nn*10-width,nn*10-width,nn*10,nn*10+width,nn*10+width),col=gene_col,border='NA')
        } else {
            polygon(c(a1[i,4],a1[i,4],a1[i,4]-0.85*(a1[i,4]-a1[i,3]),a1[i,3],a1[i,4]-0.85*(a1[i,4]-a1[i,3]),a1[i,4]), c(nn*10+width,nn*10-width,nn*10-width,nn*10,nn*10+width,nn*10+width),col=gene_col,border='NA')
        }
    }
}

axis(1)
