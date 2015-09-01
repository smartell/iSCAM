#===============================================================
#Maps to illustrate the TMA
#
#===============================================================


##if(!require('rgdal'))         install.packages('rgdal')
##if(!require('mapplots'))         install.packages('mapplots')
##
##
##library("rgdal")
##library(ggplot2)
##library(ggmap)
##library(scales)
##library(grid)
##library(reshape)
##library(maptools)
##library(mapplots)
##
##setwd("/Users/catarinawor/Documents/iSCAM/src/mseapp/data/hp/shapefiles/Reg_Area_shapefiles")
##
##US <-readOGR("/Users/catarinawor/Documents/iSCAM/src/mseapp/data/hp/shapefiles/Reg_Area_shapefiles","IPHC_RegulatoryAreas_US")
##US <- spTransform(US, CRS("+proj=longlat +datum=WGS84"))
##US <- fortify(US)
##
##CAN <-readOGR(".","IPHC_RegulatoryAreas_CAN")
##CAN <- spTransform(CAN, CRS("+proj=longlat +datum=WGS84"))
##CAN <- fortify(CAN)
##
##
##lat<-rep(c(42.5,51,55.5,58,54,52,50,62),2)
##long<-rep(c(-126.5,-131,-137,-145,-157,-167,-175,-168),2)
##value<-c(c(0.87,5.77,5.48,12.58,4.04,2.04,1.29,4.34),c(1.07,7.73,5.99,13.59,4.54,1.87,1.57,6.17))
##labels<-rep(c("2A","2B","2C","3A","3B","4A","4B","4CDE"),2)
##type<-c(rep("Blue line",8),rep("Catches",8))
##
##
##df<-data.frame(Latitude=lat, Longitude=long,value=value,type=type, area=labels)
##
##pointLabels<-annotate("text",x=df$Longitude,y=df$Latitude,size=5,font=3,fontface="bold",family="Helvetica",label=as.vector(df$area),color="white")
##
##
##basemap<-get_stamenmap(bbox = c(left = -180, bottom = 37, right = -115, top = 68),zoom = 5, maptype ="toner")
##myMap <- ggmap(basemap,base_layer=ggplot(aes(x=long,y=lat), data=US), extent =  "normal", ylab = "Latitude", xlab = "Longitude", maprange=FALSE) 
##
##
##myMap <- myMap + geom_polygon(aes(x=long, y=lat, group=group), fill='grey', size=0.5 ,color='white', data=US, alpha=0)
##myMap <- myMap + geom_polygon(aes(x=long, y=lat, group=group), fill='grey', size=0.5 ,color='white', data=CAN, alpha=0)
##myMap <- myMap + coord_map(projection="mercator", 
##              xlim=c(attr(basemap, "bb")$ll.lon, attr(basemap, "bb")$ur.lon),
##              ylim=c(attr(basemap, "bb")$ll.lat, attr(basemap, "bb")$ur.lat))
##
## 
##p2 <- myMap + geom_point(alpha=0.8,aes(x=Longitude, y=Latitude,size=value, shape=type, color=type),data=df) + labs(x = 'Longitude', y = 'Latitude')+pointLabels
##p2 <- p2 + scale_shape_manual(values=c(16,21)) + scale_fill_discrete(na.value=NA, guide="none")
##p2 <- p2 + scale_color_manual(values=c("deepskyblue1","gray80")) + scale_size_area(guide = "none", max_size = 25)
##p2 <- p2 + theme(legend.justification=c(1,1), legend.position=c(1,1)
##	,legend.background = element_rect(fill=alpha('gray80', 1))
##	,plot.background = element_rect(fill = 'white', colour = 'white'))
##p2 
##
##
##
###=========================================================================================
### pie chart attempts
##
##if(!require('maptools'))         install.packages('maptools')
##
##library(maptools)
##
###setwd("/Users/catarinawor/Documents/iSCAM/src/mseapp/data/hp/shapefiles/Reg_Area_shapefiles")
##
##US2 <-readOGR("/Users/catarinawor/Documents/iSCAM/src/mseapp/data/hp/shapefiles/Reg_Area_shapefiles","IPHC_RegulatoryAreas_US")
##US2 <- spTransform(US2, CRS("+proj=longlat +datum=WGS84"))
##
##CAN2 <-readOGR(".","IPHC_RegulatoryAreas_CAN")
##CAN2 <- spTransform(CAN2, CRS("+proj=longlat +datum=WGS84"))
##
##country<-readOGR(".","country")
##country <- spTransform(country, CRS("+proj=longlat +datum=WGS84"))
##
##setwd("/Users/catarinawor/Documents/iSCAM/src/mseapp/data/hp")
##catdat<-read.csv("2014_catch.csv")
##
##
##lat<-rep(c(42.5,51,55.5,58,54,52,50,62),3)
##long<-rep(c(-126.5,-131,-137,-145,-157,-167,-175,-168),3)
##
##df2<-cbind(catdat[catdat$Type=="Catches",],lat,long)
##
##df3<-df2[df2$Type=="Catches",]
##
##
##
##lat1<-c(42.5,51,55.5,58,54,52,50,62)
##long1<-c(-126.5,-131,-137,-145,-157,-167,-175,-168)
##labels1<-c("2A","2B","2C","3A","3B","4A","4B","4CDE")
##
##
##plot(US2,xlim=c(-180,-115),ylim=c(37,65))  
##plot(CAN2,xlim=c(-180,-115),ylim=c(37,65), add=T)
##plot(country,xlim=c(-180,-115),ylim=c(37,65),add=T, col="grey40") 
##
##box()
##axis(1, at = seq(-180,-115,by=5), labels = TRUE, tick = TRUE)
##axis(2, at = seq(40,65,by=5), labels = TRUE, tick = TRUE)
##mtext("Latitude", side = 2, line = 2.5)
##mtext("Longitude", side = 1, line = 2.5)
##text(y=lat1, x=long1 , labels = labels1, pos=1,offset = 1,cex=1.4, font=2)
##legend("bottomleft" , legend=unique(df3$sector), fill =c("olivedrab3", "deepskyblue3","darkorange2"))
##
##
##for(i in 1:8){
##  tmp<-unique(df3$Area)[i]
##  df4<-df3[df3$Area==tmp,]
##  add.pie(df4$Catch, x = df4$long, y = df4$lat , labels = " ", radius = 1, col=c("olivedrab3", "deepskyblue3","darkorange2"))
##
##}
##
##
##sum_cat<-aggregate(df3$Catch, by=list(df3$sector) , FUN=sum)
##add.pie(sum_cat$x, x = -150, y = 40 , labels = " ", radius = 2, col=c("olivedrab3", "deepskyblue3","darkorange2"))
##text(y=40, x= -150 , labels = "Total", pos=1,offset = 2,cex=1.4, font=2)








