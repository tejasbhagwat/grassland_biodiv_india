library(terra)
library(sf)
library(dplyr)
library(lubridate)
library(tidyr)
#-------------------------------------------------------------------------------------------------------------------------------------
# setwd("C:/Users/.../Desktop/.../NDVI_plotting")

#-------------------------------------------------------------------------------------------------------------------------------------
# Load the data 
#-------------------------------------------------------------------------------------------------------------------------------------
ts1 = read.csv("Chapter2/INDVI_plotting/original_ndvi_ts.csv") %>% mutate(date = mdy(system.time_start)) %>% mutate(type = "original_ndvi") %>% 
                                           dplyr::select(date, NDVI, type)

ts2 = read.csv("Chapter2/INDVI_plotting/gap_filled_ndvi_ts.csv") %>% mutate(date = mdy(system.time_start)) %>% mutate(type = "gap_filled_ndvi") %>%
                                           dplyr::select(date, NDVI, type)

ts3 = read.csv("Chapter2/INDVI_plotting/smoothed_ndvi_ts.csv") %>% mutate(date = mdy(system.time_start)) %>% mutate(type = "smoothed_ndvi") %>% 
                                           dplyr::select(date, NDVI, type) #%>% rename(NDVI = "NDVI_sg")

ndvi_ts = rbind(ts1,ts2,ts3)
#-------------------------------------------------------------------------------------------------------------------------------------
# Add cumulative NDVI value per type category
indvi_ts = ndvi_ts %>% replace(is.na(.), 0) %>% arrange(type) %>% 
           group_by(type) %>%
           mutate(indvi = cumsum(NDVI))

indvi_ts = indvi_ts %>% filter(type == "smoothed_ndvi")

library(viridis)
indvi_ts %>%
ggplot(aes(x = date, y = indvi, fill = indvi)) +
  geom_col() + 
  scale_fill_viridis_c(direction = -1) + 
  theme(legend.position =c(.9,.85), legend.title = element_blank(),
        axis.text = element_text(size = 11), legend.text = element_text(size = 11)) +
  #scale_fill_discrete(values = viridis_pal()) +
  theme_bw() +
  labs(x="", y = "INDVI") +
  theme(panel.background = element_blank(), panel.border = element_blank(),
        legend.position = "none", axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 
  #scale_x_date(date_breaks = "2 months", date_labels = "%b %Y") 

#ggsave("INDVI_smoothed.png", last_plot(), height = 1, width = 10, dpi = 400)

#-------------------------------------------------------------------------------------------------------------------------------------
# Convert 'type' to a factor with desired order
ndvi_ts$type_ro <- factor(ndvi_ts$type, levels = c("original_ndvi", "gap_filled_ndvi", "smoothed_ndvi"))

library(ggplot2)
head(ndvi_ts)
ggplot(ndvi_ts, aes(x = date, y = NDVI, color = type_ro)) + theme_bw() + 
       geom_point(size = 1.5, aes(alpha = type_ro)) + 
       geom_line(aes(alpha = type, linetype = type_ro)) + 
       theme(legend.position =c(.9,.85), legend.title = element_blank(),
             axis.text = element_text(size = 11), legend.text = element_text(size = 11)) +
       scale_color_manual(values = c("original_ndvi" = "#006400", "gap_filled_ndvi" = "#8B0000", "smoothed_ndvi" = "#4B0082")) +
       #scale_color_manual(values = c("gap_filled_ndvi" = "red", "original_ndvi" = "blue", "smoothed_ndvi" = "green")) +
       scale_alpha_manual(values = c("gap_filled_ndvi" = 0.7, "original_ndvi" = 0.3, "smoothed_ndvi" = 1)) +
       scale_linetype_manual(values = c("gap_filled_ndvi" = "dashed", "original_ndvi" = "dashed", "smoothed_ndvi" = "solid")) +
       scale_x_date(date_breaks = "2 months", date_labels = "%b %Y") +
       labs(x = "", y = "Median buffer NDVI") +
       theme(panel.background = element_rect(fill = "white"), panel.border = element_blank())
ggsave("NDVI_processing_year_round_cultivation.png", last_plot(), height = 3, width = 10, dpi = 400)
