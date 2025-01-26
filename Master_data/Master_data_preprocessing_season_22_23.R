#=========================================================================================================================================
# Bird data processing for the 2022-2023 field season
#=========================================================================================================================================
library(raster)
library(dplyr)
library(sf)
library(tidyr)
library(terra)
library(ggplot2)
library(elevatr)
library(scales)
library(stats)
library(tmap)
library(stringr)
library(vegan)
#library(MASS)
#=========================================================================================================================================
# Set the working directory
setwd("working_directory")
#=========================================================================================================================================
#-----------------------------------------------------------------------------------------------------------------------------------------
# Loading grid data 
grid_data = read_sf("Grids_2022_2023_Final_150mBuff.shp") %>% st_drop_geometry()
length(unique(grid_data$grid_name))
length(unique(grid_data$grid_id))
table(grid_data$hab_main, grid_data$grid_name)
#-----------------------------------------------------------------------------------------------------------------------------------------
# Loading bird data 
bird_data = read.csv("bird_data_master.csv") %>% dplyr::select(-c(X.1, X))
head(bird_data)

#=========================================================================================================================================
# Loading mean elevation and mean NDVI information
#=========================================================================================================================================
#' Elevation can be alternatively also be done using the 'elevatr' package
grid_elev = read.csv("grids_mean_elevation.csv") %>% dplyr::select(grid_id, mean_elevation)

#=========================================================================================================================================
# Working with all species 
#=========================================================================================================================================
#-----------------------------------------------------------------------------------------------------------------------------------------
# Explore species distribution
unique(sort(bird_data$species))
length(unique(bird_data$grid_id))
table(bird_data$hab_main)
table(bird_data$grid_id)
table(bird_data$species)

# Number of transects per species
bird_data_n = bird_data %>% group_by(species) %>% summarise(Abundance = sum(n))
#-----------------------------------------------------------------------------------------------------------------------------------------
# Importing species names to have a consistent names column
#obc_df = read.csv("bird_data_names_final.csv")
#obc_df = obc_df %>% mutate(species = Common.name)
#obc_df$species = gsub("-"," ", obc_df$species)
#obc_df = left_join(obc_df, bird_data_n)

#write.csv(obc_df, "bird_data_names_final_2.csv")
# Generating species names
species_names = unique(c(bird_data$species))
#-----------------------------------------------------------------------------------------------------------------------------------------
# Rearranging the data.frame
head(bird_data)
bird_data_form = bird_data %>% group_by(species, grid_id) %>% 
                               summarise(n_total = sum(n)) %>%
                               pivot_wider(id_cols = "grid_id", names_from = "species", values_from = "n_total") %>%
                               dplyr::select(grid_id, all_of(species_names))
# Changing NAs to 0 
bird_data_form[is.na(bird_data_form)] = 0

#-----------------------------------------------------------------------------------------------------------------------------------------
# Adding overall abundance and richness
bird_data_form$overall_abundance = rowSums(bird_data_form[2:145])
bird_data_form$overall_richness = specnumber(bird_data_form[2:145])
#write.csv(bird_data_form, "C:/Tejas_Uni_Goe/PhD/04_Data/02_India/05_Modeling/overall_abun_rich.csv")
#=========================================================================================================================================
# Adding Migration status and habitat guild (for migratory birds only)
#=========================================================================================================================================
# Separate abundance-richness responses based on migratory and non-migratory birds (Could be modeled separately)
LD_migrant_names = c("Yellow wagtail","Short toed lark","Tawny pipit","Citrine wagtail","Desert wheatear","Eurasian Hoopoe",
                     "Isabelline wheatear","Blyth's reed warbler","Booted warbler","Siberian stonechat","Rosy starling",
                     "Bluethroat","Steppe eagle","Common kestrel","Short eared owl","Black headed bunting","Red headed bunting",
                     "Montegu's harrier","Pallid harrier","Grey necked bunting","Marsh harrier","Common quail", "Common rosefinch",
                     "Eastern imperial eagle", "Eurasian sparrowhawk", "Grey wagtail")

length(unique(LD_migrant_names)) # 26 migratory species
bird_data$move_status = NA
bird_data$move_status = ifelse(bird_data$species %in% LD_migrant_names, "Long-distance migrant","Asia-subcontinent migrant")
table(bird_data$move_status)
#unique(bird_data$grid_id,bird_data$move_status)
#write.csv(bird_data, "C:/Tejas_Uni_Goe/PhD/04_Data/02_India/02_Field_Data/01_Master_Birds/Season_2022_23/bird_data_master_updated2.csv")
LD_migrant_names_df = data.frame(species = LD_migrant_names)
#write.csv(LD_migrant_names_df, "C:/Tejas_Uni_Goe/PhD/04_Data/02_India/02_Field_Data/01_Master_Birds/Season_2022_23/LD_species_names.csv")

#-----------------------------------------------------------------------------------------------------------------------------------------
# Seprate the data frame based on the movement status
LD_migrants = bird_data %>% dplyr::filter(move_status == "Long-distance migrant")
#-----------------------------------------------------------------------------------------------------------------------------------------
# Formatting long-distance migrant data frame for modeling
LD_migrants_model_df = LD_migrants %>% group_by(species, grid_id) %>% 
                                       summarise(n_total = sum(n)) %>%
                                       pivot_wider(id_cols = "grid_id", names_from = "species", values_from = "n_total")
LD_migrants_model_df[is.na(LD_migrants_model_df)] = 0
colnames(LD_migrants_model_df)
#sum(LD_migrants_model_df[which(LD_migrants_model_df$grid_id == "G110.7"),c(2:27)])
#-----------------------------------------------------------------------------------------------------------------------------------------
# Adding long-distance migrants richness information
LD_migrants_model_df$LD_richness = specnumber(LD_migrants_model_df[c(2:27)])
LD_migrants_model_df = LD_migrants_model_df %>% dplyr::select(grid_id, LD_richness)
#plot(density(LD_migrants_model_df$richness))
#-----------------------------------------------------------------------------------------------------------------------------------------
# Formatting subcontinent birds data frame for modeling
subcontinent = bird_data %>% dplyr::filter(move_status == "Asia-subcontinent migrant")
subcontinent_model_df = subcontinent %>% group_by(species, grid_id) %>% 
                                         summarise(n_total = sum(n)) %>%
                                         pivot_wider(id_cols = "grid_id", names_from = "species", values_from = "n_total")
subcontinent_model_df[is.na(subcontinent_model_df)] = 0
colnames(subcontinent_model_df)
#sum(subcontinent_model_df[which(subcontinent_model_df$grid_id == 67),c(2:119)])
# Adding richness information
subcontinent_model_df$sub_richness = specnumber(subcontinent_model_df[c(2:119)])
subcontinent_model_df = subcontinent_model_df %>% dplyr::select(grid_id, sub_richness)
#-----------------------------------------------------------------------------------------------------------------------------------------
# Separate LD_migrants further by habitat usage, load the pre-exported habitat reference (For long-distance migrants only)
#head(LD_migrants)
habitat_ref = read.csv2("long_distance_species_names.csv", sep = ",")[-1]

# Add habitat information 
LD_migrants = left_join(LD_migrants, habitat_ref)
table(LD_migrants$habitat)
#=========================================================================================================================================
# Create 5 different responses - Subcontinent, migratory, migratory - open-ground, migratory - shrub and migratory - carnivores
# Output will be saved in "Abundance richness models" folder   
#=========================================================================================================================================
#-------------------
# Subcontinent birds
subcontinent_DF = subcontinent %>% group_by(grid_id) %>% summarise(abundance = sum(n)) %>% inner_join(subcontinent_model_df)
write.csv(subcontinent_DF, "subcontinent.csv")

#----------------
# Migratory birds
LD_migrants_DF = LD_migrants %>% group_by(grid_id) %>% summarise(abundance = sum(n)) %>% inner_join(LD_migrants_model_df)
write.csv(LD_migrants_DF, "LD_migrants.csv")

#----------------
# Habitat guilds
# Migratory birds - open-ground 
LD_migrants_OG = LD_migrants %>% filter(habitat == "open-ground") %>% group_by(grid_id) %>% summarise(abundance = sum(n)) 
write.csv(LD_migrants_OG, "LD_migrants_OG.csv")

# Migratory birds - shrub
LD_migrants_SHR = LD_migrants %>% filter(habitat == "shrub") %>% group_by(grid_id) %>% summarise(abundance = sum(n)) 
write.csv(LD_migrants_SHR, "LD_migrants_SHR.csv")

# Migratory birds - carnivores
LD_migrants_CAR = LD_migrants %>% filter(habitat == "carnivore") %>% group_by(grid_id) %>% summarise(abundance = sum(n)) 
write.csv(LD_migrants_CAR, "LD_migrants_CAR.csv")
