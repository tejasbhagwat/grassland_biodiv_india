#==================================================================================================================
# Load the necessary packages
library(terra); library(raster); library(ggplot2); library(sf); library(imager); library(tmap); library(dplyr)
library(exactextractr); library("MASS"); library(performance); library(forcats); library(bayesplot);
library(landscapemetrics); library(tidyr)
#==================================================================================================================
# Load the bird data (Files exported from Master_data_preprocessing_season_22_23.R)
#==================================================================================================================
setwd("path_to_the_wd")
grid_150m_land_use = read.csv("grid_150m_land_use.csv")
#-------------------------------
# Subcontinent/Overall response
birds_all = read.csv("overall_abun_rich.csv")[-1] %>% 
                     dplyr::select(grid_id, overall_abundance, overall_richness)
birds_sub = read.csv("subcontinent.csv")[-1]
#--------------------
# Migratory response
birds_mig = read.csv("LD_migrants.csv")[-1]
birds_mig_og = read.csv("LD_migrants_OG.csv")[-1]
birds_mig_shr = read.csv("LD_migrants_SHR.csv")[-1]
birds_mig_car = read.csv("LD_migrants_CAR.csv")[-1]

#==================================================================================================================
# Adding all the five response columns to land-use predictors dataframe
#==================================================================================================================
grid_150m_model_df = grid_150m_land_use %>% 
                     left_join(birds_all) %>% 
                     left_join(birds_sub) %>% 
                     left_join(birds_mig) %>%  
                     left_join(birds_mig_og) %>%
                     left_join(birds_mig_shr) %>%
                     left_join(birds_mig_car)
head(grid_150m_model_df)
# Replacing NAs with zero
grid_150m_model_df[is.na(grid_150m_model_df)] = 0

#===================================================================================================================
# Modeling 
#===================================================================================================================
#-------------------------------------------------------------------------------------------------------------------
# Scaling covariates before modeling  
model_df_scaled = grid_150m_model_df

# Columns to scale 
colnames(model_df_scaled)
columns_to_scale = names(model_df_scaled[,c(4:13)])

# Scale and overwrite 
model_df_scaled[columns_to_scale] <- scale(model_df_scaled[columns_to_scale])

#-------------------------------------------------------------------------------------------------------------------
# Check for correlation
library(ggcorrplot)
#colnames(model_df_scaled)
model_df_scaled_corr = model_df_scaled[,c(4:13)] %>% rename(semi_perennial_area_ha = intensive_area_ha, 
                                                            annual_area_ha = seasonal_area_ha,
                                                            woody_vegetation_area_ha = swf_area_ha, 
                                                            savanna_area_ha = open_area_ha,
                                                            shannon_index = shdi, mean_field_size_m2 = mean_field_size_m2,
                                                            built_up_m2 = built_up_m2, mean_elevation_m = mean_elevation_m, 
                                                            total_edge_m = total_edge_m, edge_density = edge_density)
ggcorrplot(cor(model_df_scaled_corr, method = "pearson"), method = "square", hc.order = TRUE, outline.col = "white", 
           ggtheme = ggplot2::theme_minimal, lab = TRUE, lab_size = 3, type = "full", colors = c("#6D9EC1", "white", "#E46726"),
           insig = "blank")

#ggsave("C:/Tejas_Uni_Goe/PhD/08_Writing/02_Chapter2/Final_submission/Revision/Resubmission/Figures/Appendix/LCC_pear_corr.png", last_plot(),
#       dpi = 300, height = 7, width = 7)
#-------------------------------------------------------------------------------------------------------------------
# Modeling framework
#-------------------------------------------------------------------------------------------------------------------
library(lme4)
library(sjstats)
library(sjPlot)
library(glmmTMB)
library(effects)
library(rstanarm)
library(jtools)
library(MuMIn)
library(performance)
library(bayesplot)
library(modelsummary)
library(ggridges)

#------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Removing outliars (Rare congregations of birds like Skylark >300 are considered outliers)
model_df_scaled = model_df_scaled[-which(model_df_scaled$sub_abundance > 105),]
model_df_scaled = model_df_scaled[-which(model_df_scaled$LD_abundance > 200),]
#------------------------------------------------------------------------------------------------------------------------------------------------------------------
head(model_df_scaled)
grid_observers = grids_150m %>% dplyr::select(grid_id, observers)
model_df_scaled = model_df_scaled %>% left_join(grid_observers)
table(model_df_scaled$observers)
length(unique(model_df_scaled$observers))
colnames(model_df_scaled)
attach(model_df_scaled)
#===========================================================================================================================
# Overall richness and abundance models
#===========================================================================================================================
#------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Stan_glmer GLMM (All species)
#------------------------------------------------------------------------------------------------------------------------------------------------------------------
m1_stan <- stan_glmer(overall_abundance ~ intensive_area_ha + swf_area_ha + seasonal_area_ha + shdi +
                      mean_field_size_m2 + built_up_m2 + mean_elevation_m + (1|grid_name),
                      family = poisson(link = "log"),
                      data = model_df_scaled,
                      chains = 4,     # Number of Markov chains
                      cores = 3,      # Number of CPU cores to use
                      iter = 6000,    # Number of iterations
                      warmup = 2000,  # Number of warmup iterations
                      seed = 123      # Set a random seed for reproducibility
)
m2_stan <- stan_glmer(overall_richness ~ intensive_area_ha + swf_area_ha + seasonal_area_ha + shdi +
                      mean_field_size_m2 + built_up_m2 + mean_elevation_m + (1|grid_name),
                      family = poisson(link = "log"),
                      data = model_df_scaled,
                      chains = 4,     # Number of Markov chains
                      cores = 3,      # Number of CPU cores to use
                      iter = 6000,    # Number of iterations
                      warmup = 2000,  # Number of warmup iterations
                      seed = 123      # Set a random seed for reproducibility
)
summary(m1_stan, pars = c("intensive_area_ha","swf_area_ha","seasonal_area_ha","shdi","mean_field_size_m2"), digits = 3, probs = c(0.025, 0.5, 0.975))
summary(m2_stan, pars = c("intensive_area_ha","swf_area_ha","seasonal_area_ha","shdi","mean_field_size_m2"), digits = 3, probs = c(0.025, 0.5, 0.975))
#------------------------------------------------------------------------------------------------------------------------------
# Plotting results - Creating intervals 
#------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------
# All species
#------------------------------------------------------------------------------------------------------------------------------
intervals_m1 = mcmc_intervals(m1_stan, 
                              pars = c("intensive_area_ha","swf_area_ha","mean_field_size_m2",
                                       "seasonal_area_ha","shdi"), 
                              point_est = "mean", prob = 0.50, prob_outer = 0.95, 
                              outer_size = 1.5, inner_size = 3)

intervals_m2 = mcmc_intervals(m2_stan, 
                              pars = c("intensive_area_ha","swf_area_ha","mean_field_size_m2",
                                       "seasonal_area_ha","shdi"), 
                              point_est = "mean", prob = 0.50, prob_outer = 0.95, 
                              outer_size = 1.5, inner_size = 3)
#------------------------------------------------------------------------------------------------------------------------------
# create a data.frame for all species  
#------------------------------------------------------------------------------------------------------------------------------
intervals_abun_all = data.frame(parameter = c("intensive_area_ha","swf_area_ha","mean_field_size_m2",
                                              "seasonal_area_ha","shdi"), 
                                low_end = intervals_m1$data$ll, high_end = intervals_m1$data$hh,
                                lower_bounds = intervals_m1$data$l, upper_bounds = intervals_m1$data$h, 
                                estimate = intervals_m1$data$m, significant = intervals_m1$data$l > 0 | intervals_m1$data$h < 0, 
                                response = "Abundance",
                                species_grp = "All species")
intervals_rich_all = data.frame(parameter = c("intensive_area_ha","swf_area_ha","mean_field_size_m2",
                                              "seasonal_area_ha","shdi"),
                                low_end = intervals_m2$data$ll, high_end = intervals_m2$data$hh,
                                lower_bounds = intervals_m2$data$l, upper_bounds = intervals_m2$data$h, 
                                estimate = intervals_m2$data$m, significant = intervals_m2$data$l > 0 | intervals_m2$data$h < 0, 
                                response = "Richness",
                                species_grp = "All species")

intervals_all = rbind(intervals_abun_all, intervals_rich_all)

#===========================================================================================================================
# Resident richness and abundance models
#===========================================================================================================================
#------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Stan_glmer GLMM (Species of the subcontinent)
#------------------------------------------------------------------------------------------------------------------------------------------------------------------
m3_stan <- stan_glmer(sub_abundance ~ intensive_area_ha + swf_area_ha + seasonal_area_ha + shdi +
                      mean_field_size_m2 + built_up_m2 + mean_elevation_m + (1|grid_name),
                      family = poisson(link = "log"),
                      data = model_df_scaled,
                      chains = 4,     # Number of Markov chains
                      cores = 3,      # Number of CPU cores to use
                      iter = 6000,    # Number of iterations
                      warmup = 2000,  # Number of warmup iterations
                      seed = 123      # Set a random seed for reproducibility
)
m4_stan <- stan_glmer(sub_richness ~ intensive_area_ha + swf_area_ha + seasonal_area_ha + shdi +
                      mean_field_size_m2 + built_up_m2 + mean_elevation_m + (1|grid_name),
                      family = poisson(link = "log"),
                      data = model_df_scaled,
                      chains = 4,     # Number of Markov chains
                      cores = 3,      # Number of CPU cores to use
                      iter = 6000,    # Number of iterations
                      warmup = 2000,  # Number of warmup iterations
                      seed = 123      # Set a random seed for reproducibility
)
summary(m3_stan, pars = c("intensive_area_ha","swf_area_ha","seasonal_area_ha","shdi","mean_field_size_m2"), digits = 3, probs = c(0.025, 0.5, 0.975))
summary(m4_stan, pars = c("intensive_area_ha","swf_area_ha","seasonal_area_ha","shdi","mean_field_size_m2"), digits = 3, probs = c(0.025, 0.5, 0.975))
#------------------------------------------------------------------------------------------------------------------------------
# Resident species
#------------------------------------------------------------------------------------------------------------------------------
intervals_m3 = mcmc_intervals(m3_stan, 
                              pars = c("intensive_area_ha","swf_area_ha","mean_field_size_m2",
                                       "seasonal_area_ha","shdi"), 
                              point_est = "mean", prob = 0.50, prob_outer = 0.95, 
                              outer_size = 1.5, inner_size = 3)

intervals_m4 = mcmc_intervals(m4_stan, 
                              pars = c("intensive_area_ha","swf_area_ha","mean_field_size_m2",
                                       "seasonal_area_ha","shdi"), 
                              point_est = "mean", prob = 0.50, prob_outer = 0.95, 
                              outer_size = 1.5, inner_size = 3)
#------------------------------------------------------------------------------------------------------------------------------
# create a data.frame for all species  
#------------------------------------------------------------------------------------------------------------------------------
intervals_abun_all = data.frame(parameter = c("intensive_area_ha","swf_area_ha","mean_field_size_m2",
                                              "seasonal_area_ha","shdi"), 
                                low_end = intervals_m3$data$ll, high_end = intervals_m3$data$hh,
                                lower_bounds = intervals_m3$data$l, upper_bounds = intervals_m3$data$h, 
                                estimate = intervals_m3$data$m, significant = intervals_m3$data$l > 0 | intervals_m3$data$h < 0, 
                                response = "Abundance",
                                species_grp = "All species")
intervals_rich_all = data.frame(parameter = c("intensive_area_ha","swf_area_ha","mean_field_size_m2",
                                              "seasonal_area_ha","shdi"),
                                low_end = intervals_m4$data$ll, high_end = intervals_m4$data$hh,
                                lower_bounds = intervals_m4$data$l, upper_bounds = intervals_m4$data$h, 
                                estimate = intervals_m4$data$m, significant = intervals_m4$data$l > 0 | intervals_m4$data$h < 0, 
                                response = "Richness",
                                species_grp = "All species")

intervals_all = rbind(intervals_abun_all, intervals_rich_all)

#====================================================================================================================================
# Migratory bird species model
#====================================================================================================================================
colnames(model_df_scaled)
plot(density(model_df_scaled$abundance_mig_og))
table(model_df_scaled$abundance_mig_car)
#-------------------------------
# Stan_glmer GLMM
LDM1_stan <- stan_glmer(LD_abundance ~ intensive_area_ha + swf_area_ha + seasonal_area_ha + shdi + mean_field_size_m2 +
                        mean_elevation_m + built_up_m2 + (1|grid_name),
                        family = poisson(link = "log"),
                        data = model_df_scaled,
                        chains = 4,     # Number of Markov chains
                        cores = 3,      # Number of CPU cores to use
                        iter = 6000,    # Number of iterations
                        warmup = 2000,  # Number of warmup iterations
                        seed = 123      # Set a random seed for reproducibility
)

LDM2_stan <- stan_glmer(LD_richness ~ intensive_area_ha + swf_area_ha + seasonal_area_ha + shdi + mean_field_size_m2 + 
                        mean_elevation_m + built_up_m2 + (1|grid_name),
                        family = poisson(link = "log"),
                        data = model_df_scaled,
                        chains = 4,     # Number of Markov chains
                        cores = 3,      # Number of CPU cores to use
                        iter = 6000,    # Number of iterations
                        warmup = 2000,  # Number of warmup iterations
                        seed = 123      # Set a random seed for reproducibility
)

summary(LDM1_stan, pars = c("intensive_area_ha","swf_area_ha","seasonal_area_ha","shdi","mean_field_size_m2","mean_elevation_m",
                            "built_up_m2"), digits = 3, probs = c(0.025, 0.5, 0.975))
summary(LDM2_stan, pars = c("intensive_area_ha","swf_area_ha","seasonal_area_ha","shdi","mean_field_size_m2","mean_elevation_m",
                            "built_up_m2"), digits = 3, probs = c(0.025, 0.5, 0.975))
#---------------------------------
# saving bayesian output
intervals_LDM1 = mcmc_intervals(LDM1_stan, pars = c("intensive_area_ha","swf_area_ha","mean_field_size_m2",
                                                    "seasonal_area_ha","shdi"),
                                point_est = "mean", prob = 0.50, prob_outer = 0.95, outer_size = 1.2, inner_size = 2.0)

intervals_LDM2 = mcmc_intervals(LDM2_stan, pars = c("intensive_area_ha","swf_area_ha","mean_field_size_m2",
                                                    "seasonal_area_ha","shdi"),
                                point_est = "mean", prob = 0.50, prob_outer = 0.95, outer_size = 1.2, inner_size = 2.0)
#---------------------------------
# create a data.frame of results 
intervals_LDM_abun = data.frame(parameter = c("intensive_area_ha","swf_area_ha","mean_field_size_m2",
                                              "seasonal_area_ha","shdi"),
                            low_end = intervals_LDM1$data$ll, high_end = intervals_LDM1$data$hh,
                            lower_bounds = intervals_LDM1$data$l, upper_bounds = intervals_LDM1$data$h, 
                            estimate = intervals_LDM1$data$m, significant = intervals_LDM1$data$l > 0 | intervals_LDM1$data$h < 0, 
                            response = "Abundance",
                            species_grp = "Long distance migratory")
intervals_LDM_rich = data.frame(parameter = c("intensive_area_ha","swf_area_ha","mean_field_size_m2",
                                              "seasonal_area_ha","shdi"),
                            low_end = intervals_LDM2$data$ll, high_end = intervals_LDM2$data$hh,
                            lower_bounds = intervals_LDM2$data$l, upper_bounds = intervals_LDM2$data$h, 
                            estimate = intervals_LDM2$data$m, significant = intervals_LDM2$data$l > 0 | intervals_LDM2$data$h < 0, 
                            response = "Richness",
                            species_grp = "Long distance migratory")
intervals_LDM_all = rbind(intervals_LDM_abun, intervals_LDM_rich)

#--------------------------------------------------
# Combining the ALL SPECIES with MIGRATORY SPECIES 
intervals_all_birds = rbind(intervals_all, intervals_LDM_all)

#====================================
# Create the custom ggplot2 plot
# Create the plot
alph = ifelse(intervals_all_birds$significant == T, 1, 0.20)
pos <- position_nudge(y = ifelse(intervals_all_birds$response == "Abundance", 0.1, -0.1))

# Plotting intervals for all birds
intervals_all_birds_plot = intervals_all_birds %>% 
                           mutate(parameter = fct_relevel(parameter, "mean_field_size_m2",
                                                          "swf_area_ha", "seasonal_area_ha",
                                                          "intensive_area_ha", "shdi")) %>%
                           mutate(species_grp = recode(species_grp, "All species" = "Birds of the Indian sub-continent \n(n=118)", "Long distance migratory" = "Palearctic migratory birds \n(n=26)"),
                                  response = recode(response, "Abundance" = "Abundance", "Richness" = "Species richness")) %>%
                           ggplot(aes(x = estimate, y = parameter, color = response, alpha = alph)) + 
                           scale_color_manual(values = c('#800000', '#3CB371')) + scale_alpha_identity() + 
                           scale_y_discrete(labels = c("shdi" = "Land cover diversity",
                                                       "intensive_area_ha" = "Proportion of semi-perennial crops",
                                                       "seasonal_area_ha" = "Proportion of annual crops",
                                                       "swf_area_ha" = "Proportion of woody vegetation",
                                                       "mean_field_size_m2" = "Mean field size")) +
                           geom_linerange(aes(xmin = lower_bounds, xmax = upper_bounds), position = pos, linewidth = 2) +
                           geom_linerange(aes(xmin = low_end, xmax = high_end), position = pos) + theme_minimal() +
                           geom_point(position = pos, color="black", size = 2.1) + geom_vline(xintercept=0, linetype="dashed") + 
                           labs(title = "") + facet_grid(~species_grp) +
                           theme(legend.position = "bottom", legend.title = element_blank(), 
                                 plot.background = element_rect(fill = 'white', colour = "white"),
                                 strip.text = element_text(size = 16), legend.text = element_text(size = 16),
                                 axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 16),
                                 axis.title.x = element_blank(), axis.title.y = element_blank(),
                                 plot.margin = margin(10, 10, 10, 10))
ggsave("Overall_models_26_11_24.png", intervals_all_birds_plot, dpi = 300, height = 5, width = 10.5)

#=======================================================
# Guild Models
#=======================================================
#' Open area / grassland class is strongly correlated with 5 different variables. So just replacing
#' intensive area with open does not work because then open area and seasonal area are approx. 60 correlated and 
#' that brings has confounding effects in the model?
#' 
#g1 = glmmTMB(abundance_mig_car ~ intensive_area_ha + mean_field_size_m2 + swf_area_ha + 
#             seasonal_area_ha + shdi + mean_field_size_m2:swf_area_ha + 
#             (1 | grid_id), family = poisson, data = model_df_scaled)
#summary(g1)
#g2 = glmmTMB(abundance_mig_og ~ intensive_area_ha + mean_field_size_m2 + swf_area_ha + 
#             seasonal_area_ha + shdi + mean_field_size_m2:swf_area_ha + 
#               (1 | grid_id), family = poisson, data = model_df_scaled)
#summary(g2)
#g3 = glmmTMB(abundance_mig_shr ~ intensive_area_ha + mean_field_size_m2 + swf_area_ha + 
#             seasonal_area_ha + shdi + mean_field_size_m2:swf_area_ha + 
#               (1 | grid_id), family = poisson, data = model_df_scaled)
#summary(g3)

# Models with grassland habitats
#g1_grass = glmmTMB(abundance_mig_car ~ open_area_ha + swf_area_ha + mean_field_size_m2 + mean_elevation_m + built_up_m2 + 
#                     (1 | grid_id), family = poisson, data = model_df_scaled)
#summary(g1_grass)

#g2_grass = glmmTMB(abundance_mig_og ~ open_area_ha + swf_area_ha + mean_field_size_m2 + mean_elevation_m + built_up_m2 + 
#                  (1 | grid_id), family = poisson, data = model_df_scaled)
#summary(g2_grass)

#g3_grass = glmmTMB(abundance_mig_shr ~ open_area_ha + swf_area_ha + mean_field_size_m2 + mean_elevation_m + built_up_m2 + 
#                  (1 | grid_id), family = poisson, data = model_df_scaled)
#summary(g3_grass)

#----------------------------------------------------------------
# Frequentist approach, model checking corner..
#model_diag = simulateResiduals(fittedModel = g1, plot = T)
#dispersion_test = DHARMa::testDispersion(model_diag)
#check_overdispersion(m1)

# prevent fitting sub-models to different datasets
#oop <- options(na.action = "na.fail")
#dd_3 <- dredge(g1)

#library(rstanarm); library(bayesplot)
#prop_zero <- function(y) mean(y == 0)
#pp_check(g3_stan, plotfun = "stat", stat = "prop_zero", binwidth = .005)
#pp_check(m1_stan, plotfun = "scatter_avg")
#pp_check(m1_stan, plotfun = "intervals")

#----------------------------------------------------------------
#-------------------------------
# Stan_glmer GLMM
g1_stan <- stan_glmer(abundance_mig_car ~ intensive_area_ha + mean_field_size_m2 + swf_area_ha + 
                      seasonal_area_ha + shdi + mean_elevation_m + built_up_m2 + (1|grid_name),
                      family = poisson(link = "log"),
                      data = model_df_scaled,
                      chains = 4,     # Number of Markov chains
                      cores = 3,      # Number of CPU cores to use
                      iter = 6000,    # Number of iterations
                      warmup = 2000,  # Number of warmup iterations
                      seed = 123      # Set a random seed for reproducibility
)
summary(g1_stan, pars = c("intensive_area_ha","swf_area_ha","seasonal_area_ha","shdi","mean_field_size_m2"), digits = 3, probs = c(0.025, 0.5, 0.975))
g2_stan <- stan_glmer(abundance_mig_og ~ intensive_area_ha + mean_field_size_m2 + swf_area_ha + 
                      seasonal_area_ha + shdi + mean_elevation_m + built_up_m2 + (1|grid_name),
                      family = poisson(link = "log"),
                      data = model_df_scaled,
                      chains = 4,     # Number of Markov chains
                      cores = 3,      # Number of CPU cores to use
                      iter = 6000,    # Number of iterations
                      warmup = 2000,  # Number of warmup iterations
                      seed = 123      # Set a random seed for reproducibility
)
summary(g2_stan, pars = c("intensive_area_ha","swf_area_ha","seasonal_area_ha","shdi","mean_field_size_m2"), digits = 3, probs = c(0.025, 0.5, 0.975))
g3_stan <- stan_glmer(abundance_mig_shr ~ intensive_area_ha + mean_field_size_m2 + swf_area_ha + 
                      seasonal_area_ha + shdi + mean_elevation_m + built_up_m2 + (1|grid_name),
                      family = poisson(link = "log"),
                      data = model_df_scaled,
                      chains = 4,     # Number of Markov chains
                      cores = 3,      # Number of CPU cores to use
                      iter = 6000,    # Number of iterations
                      warmup = 2000,  # Number of warmup iterations
                      seed = 123      # Set a random seed for reproducibility
)
summary(g3_stan, pars = c("intensive_area_ha","swf_area_ha","seasonal_area_ha","shdi","mean_field_size_m2"), digits = 3, probs = c(0.025, 0.5, 0.975))

#---------------------------------
# saving bayesian output
intervals_g1_stan = mcmc_intervals(g1_stan, pars = c("intensive_area_ha","swf_area_ha","seasonal_area_ha",
                                                     "shdi","mean_field_size_m2"),
                                   point_est = "mean", prob = 0.50, prob_outer = 0.95, outer_size = 1, inner_size = 2)

intervals_g2_stan = mcmc_intervals(g2_stan, pars = c("intensive_area_ha","swf_area_ha","seasonal_area_ha",
                                                     "shdi","mean_field_size_m2"), 
                                   point_est = "mean", prob = 0.50, prob_outer = 0.95, outer_size = 1, inner_size = 2)

intervals_g3_stan = mcmc_intervals(g3_stan, pars = c("intensive_area_ha","swf_area_ha","seasonal_area_ha",
                                                     "shdi","mean_field_size_m2"),
                                   point_est = "mean", prob = 0.50, prob_outer = 0.95, outer_size = 1, inner_size = 2)
#---------------------------------
# create a data.frame of results 
intervals_g1_df = data.frame(parameter = intervals_g1_stan$data$parameter, 
                             low_end = intervals_g1_stan$data$ll, high_end = intervals_g1_stan$data$hh,
                             lower_bounds = intervals_g1_stan$data$l, upper_bounds = intervals_g1_stan$data$h, 
                             estimate = intervals_g1_stan$data$m, significant = intervals_g1_stan$data$l > 0 | intervals_g1_stan$data$h < 0, 
                             species_grp = "carnivores")


intervals_g2_df = data.frame(parameter = intervals_g2_stan$data$parameter, 
                             low_end = intervals_g2_stan$data$ll, high_end = intervals_g2_stan$data$hh,
                             lower_bounds = intervals_g2_stan$data$l, upper_bounds = intervals_g2_stan$data$h, 
                             estimate = intervals_g2_stan$data$m, significant = intervals_g2_stan$data$l > 0 | intervals_g2_stan$data$h < 0, 
                             species_grp = "open-ground")

intervals_g3_df = data.frame(parameter = intervals_g3_stan$data$parameter, 
                             low_end = intervals_g3_stan$data$ll, high_end = intervals_g3_stan$data$hh,
                             lower_bounds = intervals_g3_stan$data$l, upper_bounds = intervals_g3_stan$data$h, 
                             estimate = intervals_g3_stan$data$m, significant = intervals_g3_stan$data$l > 0 | intervals_g3_stan$data$h < 0, 
                             species_grp = "shrub")

#-----------------------------------
# Combining all the guilds together
intervals_all_guilds = rbind(intervals_g1_df, intervals_g2_df, intervals_g3_df)

#----------
# Plot them
alph_guilds = ifelse(intervals_all_guilds$significant == T, 1, 0.3)
pos_guilds <- position_dodge2(width = 0.5)
#pos_guilds <- position_nudge(y = ifelse(intervals_all_birds$response == "Abundance", 0.12, -0.12))


library(forcats)
library(tidybayes)
library(RColorBrewer)
#ggplot(aes(y = parameter, x = dist_norm(estimate, )))


intervals_all_guilds %>% 
mutate(parameter = fct_relevel(parameter, "mean_field_size_m2",
                               "swf_area_ha", "seasonal_area_ha",
                               "intensive_area_ha", "shdi")) %>%
ggplot(aes(x = estimate, y = parameter, color = species_grp, alpha = alph_guilds)) + 
  scale_color_manual(values = c('#9370DB',"#8BBB7A",'#B8860B')) + 
  scale_alpha_identity() + 
  scale_y_discrete(labels = c("shdi" = "Land cover diversity",
                              "intensive_area_ha" = "Proportion of semi-perennial crops",
                              "mean_field_size_m2" = "Mean field size",
                              "seasonal_area_ha" = "Proportion of annual crops",
                              "swf_area_ha" = "Proportion of woody vegetation")) +
  geom_linerange(aes(xmin = lower_bounds, xmax = upper_bounds), position = pos_guilds, linewidth = 2) +
  geom_linerange(aes(xmin = low_end, xmax = high_end), position = pos_guilds) + theme_minimal() +
  geom_point(position = pos_guilds, color="black", size = 2.1) + geom_vline(xintercept=0, linetype="dashed") + 
  #facet_grid(~species_grp, scales = "free_x") +
  theme(legend.position = "bottom", legend.title = element_blank(), plot.title = element_text(size=16),
        axis.title.x = element_blank(), axis.title.y = element_blank(), 
        plot.background = element_rect("white", colour = "white"),
        strip.text = element_text(size = 16), legend.text = element_text(size = 16),
        axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 16)) +
  labs(title = "Abundance of the palearctic migratory birds (n=26)")

ggsave("Guild_models_26_11_2024.png", last_plot(), dpi = 300, height = 5, width = 10.5)

#=======================================================
# Posterior predictive checking
#=======================================================
y = model_df_scaled$overall_abundance # Response 
ppd = posterior_predict(m1_stan, draws = 1000) # draws from model

# Overall models
m1_stan_ppd = ppc_stat(y, ppd, stat = "mean")
#m2_stan_ppd = ppc_stat(y, ppd, stat = "mean")
#LDM1_stan_ppd = ppc_stat(y, ppd, stat = "mean")
#LDM2_stan_ppd = ppc_stat(y, ppd, stat = "mean")

# Guild models
#g1_stan_ppd = ppc_stat(y, ppd, stat = "mean")
#g2_stan_ppd = ppc_stat(y, ppd, stat = "mean")
#g3_stan_ppd = ppc_stat(y, ppd, stat = "mean")

# Plot them 
overall_models_ppd = cowplot::plot_grid(m1_stan_ppd, m2_stan_ppd, LDM1_stan_ppd, LDM2_stan_ppd, nrow = 2, ncol = 2)
guild_models_ppd = cowplot::plot_grid(g1_stan_ppd, g2_stan_ppd, g3_stan_ppd, nrow = 3, ncol = 1)
ggsave("overall_models_ppd.png", overall_models_ppd,
       dpi = 300, height = 5, width = 7)

ggsave("guild_models_ppd.png", guild_models_ppd,
       dpi = 300, height = 4, width = 5)


# Plot distribution of observed outcomes "y" and kernel density estimate of the replications of y from ppd
color_scheme_set("brightblue")
ppc_dens_overlay(y, ppd[1:50,]) + xlim(-5, 50) + theme_bw()
ppc_scatter_avg(y, ppd)
ppc_stat_2d(y, ppd, stat = c("mean","sd"))
