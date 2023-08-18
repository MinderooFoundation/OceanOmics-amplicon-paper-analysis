
library(ggplot2)
library(tidyverse)
library(mapdata)
library(sf)
library(sp)
library(rgdal)
library(terra)
library(viridis)

############
#CKI
CKI_Xlim <- c(95, 135) 
CKI_Ylim <- c(-35, -5)

CKI <- st_as_sf(map("worldHires", fill = T , plot = F, xlim = CKI_Xlim, ylim = CKI_Ylim))
Project1 <- read.csv("data/metadata/V10_CKI_eDNA_metadata_P1.csv") #just a .csv file with all our sites

png("figures/CKI_sites_JP.png", width = 15, height = 14, units = "cm", res = 300)
map_plot <- CKI %>% 
  ggplot() +
  geom_sf(fill = "grey60", lwd = 0.25)+
  coord_sf(xlim = CKI_Xlim, ylim = CKI_Ylim)+ 
  # RS
  geom_point(aes(x =118.938355 , y = -17.587892), shape = 21, fill = "#FF00FF", size = 3, stroke =0.85)+
  geom_point(aes(x = 119.318722 , y = -17.362589 ), shape = 21, fill = "#FF00FF", size = 3, stroke =0.85)+
  geom_point(aes(x =119.633606 , y =  -17.13278), shape = 21, fill = "#FF00FF", size = 3, stroke =0.85)+
  # NW WA
  geom_point(aes(x =121.963740 , y = -18.020874), shape = 24, fill = "#00FFFF", size = 3, stroke =0.85)+
  geom_point(aes(x = 122.099537 , y =-16.904656 ), shape = 24, fill = "#00FFFF", size = 3, stroke =0.85)+
  geom_point(aes(x =122.841860 , y =  -16.126522), shape = 24, fill = "#00FFFF", size = 3, stroke =0.85)+
  geom_point(aes(x =123.862021 , y = -15.137128), shape = 24, fill = "#00FFFF", size = 3, stroke =0.85)+
  geom_point(aes(x = 124.666109 , y =-14.356960 ), shape = 24, fill = "#00FFFF", size = 3, stroke =0.85)+
  geom_point(aes(x =125.409035 , y =  -13.737558), shape = 24, fill = "#00FFFF", size = 3, stroke =0.85)+
  #CKI sites
  geom_point(data = Project1, mapping = aes(x = longitude_dd, y = latitude_dd), shape = 23, fill = "#D55E00", size = 2, stroke =1.25)+ 
  theme_bw(16)+
  xlab("")+
  ylab("")+
  xlab(expression('Longitude'))+
  ylab(expression("Latitude"))
map_plot + scale_fill_viridis()
map_plot
dev.off()
