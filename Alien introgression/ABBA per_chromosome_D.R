library('pheatmap')
d <- read.table('per_chromosome_D.xls',sep='\t',header=T,row.names = 1)
par(oma=c(0,0,0,0),mar=c(0,0,0,0))
pheatmap(d,border_color = 'white',scale='none',legend=T,breaks=NA,cluster_rows=F,cluster_cols=F,cellwidth=20,cellheight=20,height=5,width=16.5,color=colorRampPalette(c("#fff3fa","#e90c59"))(10),filename = 't.pdf',show_rownames=T,angle_col=45)

#seq(0,400,40)