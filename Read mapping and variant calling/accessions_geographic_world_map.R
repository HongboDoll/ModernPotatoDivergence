
require(reshape)
require(rworldmap)
require(rworldxtra)
rice <- read.table("166_potato_geographic_coordinates.xls", head=T) #读入今天要用到的数据（请下载附件数据）


##更改颜色 potato world map 

op <- c(rgb(235,70,144,max=255),rgb(39,93,242,max=255))
#pdf('potato=world.pdf',width=8,height=7)
mapBubbles(rice,nameZSize="num",nameX="Longitude",pch=21,nameY="Latitude",mapRegion='world',nameZs =c('outgroup','wild','landrace','candolleanum'),nameZColour="Group",colourPalette=op,symbolSize=0.15,lwd=0.5,barOrient='vert',oceanCol=rgb(168,199,222,max=255),borderCol='grey',landCol="white",xlim=c(-130,-80),ylim=c(-56,66),addColourLegend=F,addLegend = FALSE) #addColourLegend = TRUE, colourLegendPos='bottomleft',legendPos = "bottomright",addLegend = FALSE,)
