library('RIdeogram')
library('rsvg')

a <- read.table('DM_v6.1_karyotype.xls', header=T, sep='\t')
b <- read.table('DM_v6.1_Altus_vernei_introgression.xls', header=T, sep='\t')
c <- read.table('DM_v6.1_NLR.xls', header=T, sep='\t')
d <- read.table('DM_v6.1_Kardent_vernei_introgression.xls', header=T, sep='\t')

### chromosome karyotype, comparasion
ideogram(karyotype = a, overlaid = b, colorset1 = c('white', "#88419d"),colorset2 = c('white', "#3acc76"),label = d,label_type = "heatmap", output = "chromosome_Altus_Kardent.svg") #label = c, 

### chromosome karyotype, single with NLR marker
ideogram(karyotype = a, overlaid = b, colorset1 = c('white', "#88419d"),colorset2 = c('white', "#8de5b9"),label = c,label_type = "marker", output = "chromosome_NLR.svg") #label = c, 

rsvg_pdf("chromosome_NLR.svg","chromosome_NLR.pdf")
rsvg_pdf("chromosome_Altus_Kardent.svg","chromosome_Altus_Kardent.pdf")
