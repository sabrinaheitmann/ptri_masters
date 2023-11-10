library(tidyverse)
library(magrittr)
library(ggmap)
library(ggrepel)
library(cowplot)

pal.region <- c("#fb832d","#016392")

garden.coords <- data.frame(Labs=c("Corvallis Garden","Hermiston Garden"),
                            Latitude=c(44.59,45.83),
                            Longitude=c(-123.2,-119.28),
                            Region=c("West","East"))

dat <- read.csv("origins.csv")

#get medians of lat lon for each population
dat.pop <- dat %>% group_by(Region,Population) %>% summarize(Latitude=median(Latitude),Longitude=median(Longitude))

#Register API - key is located in ignored file
#register_google(key = )

#make map
(map <- 
    get_map(location=c(left=-124.98,bottom=44,right=-114.29,top=50),
            zoom=7, maptype = "watercolor",source="osm",force=F) %>%
    ggmap()+
    geom_point(data=dat.pop,aes(x=Longitude,y=Latitude, fill=Region, shape=Region), size=5)+
    geom_point(data=garden.coords,inherit.aes = F,aes(x=Longitude,y=Latitude, shape=Region),size=5,fill="white",
               show.legend = F )+
    geom_label_repel(data=garden.coords,inherit.aes = F,aes(x=Longitude,y=Latitude,label=Labs),
                     point.padding = 0.5,box.padding = 0.6,size=6)+
    scale_fill_manual(values=pal.region)+
    scale_shape_manual(values=c(23,21))+
    labs(x="Longitude",y="Latitude")+
    # theme_map()+
    theme())


#add inset figure
ggdraw()+
  draw_plot(map)+
  draw_image("rangemap.png",scale=1,width=0.30,x=0.61,y=0.31)

#save to file
ggsave("map_with_legend.jpg",width=20,height=16,units="cm")
