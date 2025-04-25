rm(list = ls())

library(data.table)
library(ggplot2)
library(ggstance)
library(ggridges)
library(ggpubr)
library(MASS)
library(stringr)
library(scales)
library(units)
library(parallel)

sf::sf_use_s2(FALSE)

##########################################
# Reported Crime Trends by Policy Status #
##########################################

df <- fread("./data/processed/df_crime_analysis.csv")

# Plot average effects within non-treated sites
# and effects within treated sites. 

df_avg <- df[total_population >= 10000 & !is.na(burglary_total) & burglary_total != 0 & outlier == F, ]
df_avg[is.na(post_treat), post_treat := 0]
df_avg[, post_treat := as.factor(post_treat)]

# Remove tucson from 2021 (which didn't appear to report
# data to ucr). 
df_avg <- df_avg[!(year == 2021 & location == 'Tucson')]

df_avg[, burglary_count := (burglary_total*total_population)/1000]
df_avg[, burglary_count := sum(.SD$burglary_count), by = c("year", "post_treat")] 
df_avg[, pop_group := sum(.SD$total_population), by = c("year", "post_treat")]
df_avg[, burglary_total := 10000*(burglary_count/pop_group)]

df_avg[, assault_count := (assault_total*total_population)/1000]
df_avg[, assault_count := sum(.SD$assault_count), by = c("year", "post_treat")] 
df_avg[, pop_group := sum(.SD$total_population), by = c("year", "post_treat")]
df_avg[, assault_total := 10000*(assault_count/pop_group)]

df_avg <- unique(df_avg[, .(year, burglary_total, assault_total, post_treat)])
df_avg[, post_treat := ifelse(post_treat == 0, "Cities without a program", "Cities with Crime-Free Housing")]

df_avg <- melt(df_avg, id.vars = c("year", "post_treat"), value.name = "totals", variable.name= "crime_type")
df_avg[, crime_type := ifelse(crime_type == "burglary_total", "Burglaries", "Assaults")]

df_avg <- unique(df_avg)
  
p1 <- ggplot(df_avg, aes(x = year, y = totals, color = post_treat)) +
  geom_line(linewidth = 1) +
  facet_wrap(~crime_type, ncol = 2) +
  scale_x_continuous(limits = c(2010, 2020), breaks = seq(2010, 2020, 5)) +
  labs(title = "Assaults and burglaries do not decrease more in cities with crime-free housing \nthan cities without programs",
       x = "Year", 
       y = "Reported crimes per\n 10,000 people",
       color = "",
       caption = "Source: FBI UCR, 2010 - 2020. \nNotes: Burglaries include forced entry, non-forced entry, and attempted. \nAssaults include those with a gun, with a knife, with a weapon, \nunarmed, and simple assaults. Data displayed for cities in AZ\n w/ populations greater than 10,000.") +
  theme_bw() +
  theme(strip.background = element_blank(),
        legend.position = "bottom",
        axis.text = element_text(family = 'sans', size = 10),
        plot.title = element_text(family = 'sans', size = 14),
        axis.title = element_text(family = 'sans', size = 12),
        strip.text = element_text(size = 12, family = 'sans'))

ggsave(plot = p1, filename = "./figs/reported_crime_trends.pdf", height = 8.27, width = 8.27)

# Strip plot down to essential elements: change in crime between
# 2010 and 2023:

df_min <- df_avg[year == 2010 | year == 2020,]

ggplot(df_min, aes(x = year, y = totals, color = post_treat)) +
  geom_point(size = 4) +
  geom_line(size = 1) +
  labs(y = "") +
  lims(y = c(0, 200)) +
  scale_x_continuous(breaks = c(2010, 2020)) +
  scale_color_manual(values = c("#8da0cb", "#66c2a5")) +
  facet_wrap(~crime_type, ncol = 1, strip.position = "left") +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text = element_text(angle = 0),
        axis.ticks = element_blank(),
        legend.position = "none",
        line = element_blank(),
        axis.text.y = element_blank())

ggsave(filename = "./figs/crime_trend_proto.svg", height = 4.27, width = 4.27)

##################
# SFD ATT effects#
##################

# Using Rubin's rules, simulate coefficient values from each model's 
# posterior distribution and collapse estimate draws across map angles. 

# Then, display ATT effects as percentage of eviction rate for treated locs

# Read in cities and combine into single file:
combine_sfd_df <- function(file){
  
  d <- st_read(file)
  
  setDT(d)
  return(d)

}

sfd_datasets <- paste0("./data/processed/sfd/", list.files("./data/processed/sfd/"))
sfd_models   <- paste0("./data/models/", list.files("./data/models/"))

sfd_datasets <- sfd_datasets[!(sfd_datasets %like% "crime")]
sfd_models <- sfd_models[!(sfd_models %like% "crime")]

df_sfd <- setDT(ldply(sfd_datasets, combine_sfd_df))
df_sfd[, geometry := NULL]

cfh_cities <- c("avondale", "chandler", "mesa", "scottsdale")

# Read in a model for a given city. Using estimated betas and vcov,
# simulate 1000 parameter draws from each model and combine into one data.frame

sim_betas <- function(a, mod, city){
  
  submod <- mod[[a]]
  
  vars  <- submod[['coef']]$var
  
  k <- length(vars)
  
  ests  <- submod[['coef']]$Estimate
  vcovs <- submod[['vcov']][1:k, 1:k]
  
  sims <- as.data.table(mvrnorm(1000, ests, vcovs))
  
  # Second column is always the treatment variable
  treatname <- vars[2]
  
  if (length(vars) > 2){
    setnames(sims, names(sims), c(vars[1], "est_beta", vars[3:length(vars)]))
  }else{
    setnames(sims, names(sims), c(vars[1], "est_beta"))
    
  }
  sims[, treat_name := treatname]
  sims[, angle := a]
  
  return(sims)
  
}

combine_sims <- function(city, mod){
  
  mod <- mod[[city]]
  
  angles <- names(mod)[names(mod) %like% 'angle']
  
  sims <- rbindlist(lapply(angles, sim_betas, mod = mod, city = city))
  
  sims[, angle := gsub("angle=", "", angle)]
  sims[, angle := factor(angle, levels = unique(angle))]
  
  sims[, location := city]
  
  return(sims)
  
}

sim_wrapper <- function(model){
  
  mod <- readRDS(paste0(model))
  cities <- names(mod)
  
  sim_cities <- rbindlist(lapply(cities, combine_sims, mod = mod))
  
  # Make model names human-readable
  model_name <- gsub("./data/models/|.Rds|.rds", "", model)
  model_name <- gsub("_", " ", model_name)
  
  model_name <- ifelse(model_name %like% "any", 
                       gsub("cfh any", "CFH any,", model_name),
                       gsub("cfh count", "CFH count,", model_name))
  
  model_name <- ifelse(model_name %like% "adj", 
                       gsub("adj", "adjusted", model_name),
                       gsub("unadj", "unadjusted", model_name))
  
  sim_cities[, model_name := model_name]
  
  return(sim_cities)
  
}

sim_results <- rbindlist(lapply(sfd_models, sim_wrapper), fill = T)

sim_att <- sim_results[, .(est_beta, treat_name, location, model_name)]

# Calculate lower ui, mean, and upper ui estimates for each model (adjust/unadjust),
# treatment-parameter (any v. count), and city:

sim_att[, est_mean  := quantile(.SD$est_beta, probs = 0.5), by = c("model_name", "location")]
sim_att[, est_lower := quantile(.SD$est_beta, probs = 0.025), by = c("model_name", "location")]
sim_att[, est_upper := quantile(.SD$est_beta, probs = 0.975), by = c("model_name", "location")]

# Also calculate across all sites:
sim_att_overall <- sim_results[, .(est_beta, treat_name, location, model_name)]
sim_att_overall[, location := "Across Sites"]

sim_att_overall[, est_mean  := quantile(.SD$est_beta, probs = 0.5), by = c("model_name")]
sim_att_overall[, est_lower := quantile(.SD$est_beta, probs = 0.025), by = c("model_name")]
sim_att_overall[, est_upper := quantile(.SD$est_beta, probs = 0.975), by = c("model_name")]

sim_att <- unique(sim_att[, .(location, treat_name, model_name, est_mean, est_lower, est_upper)])
sim_att_overall <- unique(sim_att_overall[, .(location, treat_name, model_name, est_mean, est_lower, est_upper)])

sim_att <- rbind(sim_att, sim_att_overall)

sim_att[, location := str_to_title(location)]

# Set order of models to be more intuitive for the reader (treatment-parameter type first, then
# model type):


old_model_names <- c("CFH any, unadjusted", "CFH any, adjusted", 
                     "CFH count, unadjusted", "CFH count, adjusted")

ord_mods <- c("Unadjusted model, \nbinary treatment", "Adjusted model, \nbinary treatment", 
              "Unadjusted model, \ncontinuous treatment", "Adjusted model, \ncontinuous treatment")

sim_att[, model_name := mapvalues(model_name, old_model_names, ord_mods)]
sim_att[, model_name := factor(model_name, levels = ord_mods)]

plot_colors_bw <- c("#d9d9d9", "#969696", "#737373", "#252525")

p2 <- ggplot(sim_att, aes(x = est_mean, color = model_name, y = location)) +
      geom_errorbarh(aes(xmin = est_lower, xmax = est_upper, color = factor(model_name, levels = rev(ord_mods))), size = 2, height = 0, position = position_dodgev(height = 0.4)) +
      geom_point(aes(color = factor(model_name, levels = rev(ord_mods))), size = 3, position = position_dodgev(height = 0.4)) +
      geom_vline(xintercept = 0, linetype = 'dashed', size = 1) +
      labs(title = "Estimated ATT by Model and City",
           x = "Estimated Policy Effect",
           y = "Location",
           color = "Model") +
      theme_bw() +
      scale_color_manual(values = plot_colors_bw, limits = ord_mods) +
      theme(plot.title = element_text(hjust = 0.5, family = 'sans', size = 16),
            strip.background = element_blank(),
            legend.position = "bottom",
            legend.text = element_text(family = 'sans', size = 10),
            axis.ticks = element_line(linewidth = 1),
            axis.ticks.length = unit(5.6, "points"),
            axis.title = element_text(size = 12, family = 'sans'),
            axis.title.y = element_text(size = 12, family = 'sans', angle = 0),
            axis.text = element_text(size = 10, family = 'sans'),
            axis.text.x = element_text(size = 10, family = 'sans',
                                       margin = margin(t = 5, r = 0, b = 10, l = 0)),
            legend.title = element_text(family = 'sans', size = 14))

ggsave(plot = p2, filename = "./figs/model_results/sfd/summary_att.pdf",
       device = "pdf", height = 8.3, width = 11.7, units = "in")

# Using simulated coefficients, calculate change in predicted evictions when
# cfh policies are removed:

sim_cf <- function(modelname, loc){
  
  obs  <- df_sfd[location == loc & cfh_any == 1, ]
  sims <- sim_results[location == loc & model_name == modelname,]
  
  treat_var <- unique(sims$treat_name)
  
  truth <- mean(obs$evict_rate)
  
  # Set treatment status to either equal 1 if
  # treat variable is "cfh any" or mean across blocks,
  # if treat variable is "cfh num"
  obs <- obs[, mean(get(treat_var))]
  
  # Remove extraneous columns from sims:
  sims <- sims[, .(est_beta)]
  
  # Calculate change from truth, if estimated beta held as causal effect:
  cf_change <- as.numeric(truth - obs*sims$est_beta)
  
  df_cf <- data.table("obs" = truth, "cf" = cf_change)
  df_cf[, percent_change := (obs - cf)/obs]
  
  df_cf[, location := loc]
  df_cf[, model_name := modelname]
  
  return(df_cf)
  
}

sim_cf_overall <- function(modelname){
  
  evict_observed <- 0
  evict_cf       <- 0
  
  for (city in unique(df_sfd$location)){
    
    obs  <- df_sfd[location == city, ]
    sims <- sim_results[location == city & model_name == modelname,]
    
    # Only hold onto treated locations:
    obs <- obs[cfh_any == 1,]
    
    treat_var <- unique(sims$treat_name)
    
    truth <- sum(obs$evict_rate)
    obs <- obs[, (treat_var), with = F]
    
    # Remove extraneous columns from sims:
    sims <- sims[, .(est_beta)]
    
    # Calculate change from truth, if estimated beta held as causal effect:
    cf_change <- truth - colSums(as.matrix(obs) %*% t(as.matrix(sims)))
    
    evict_observed <- evict_observed + truth
    evict_cf       <- evict_cf + cf_change
    
  }
  
  df_cf <- data.table("obs" = evict_observed, 
                      "cf" = evict_cf)
  
  df_cf[, percent_change := (obs - cf)/obs]
  
  df_cf[, location := "Across sites"]
  df_cf[, model_name := modelname]
  
  return(df_cf)
  
}

args <- expand.grid(loc = cfh_cities,
                    modelname = unique(sim_results$model_name)[!(unique(sim_results$model_name) %like% "crime")])

sim_beta_cf <- rbindlist(mapply(sim_cf, 
                                modelname = args$modelname,
                                loc = args$loc, 
                                SIMPLIFY = F))

sim_beta_cf_overall <- rbindlist(mapply(sim_cf_overall, 
                                        modelname = args$modelname,
                                        SIMPLIFY = F))

sim_beta_cf <- rbindlist(list(sim_beta_cf, sim_beta_cf_overall))

p3 <- ggplot(sim_beta_cf, aes(x = percent_change, y = model_name, fill = after_stat(x))) +
      geom_density_ridges_gradient(scale = 1.7, alpha = 0.8) +
      scale_fill_viridis_c(name = "Percentage Change in Eviction Filing Rate", option = "C") +
      facet_wrap(~location) +
      labs(x = "Eviction Filing Rate", y = "Model") +
      theme_ridges() +
      theme(strip.background = element_blank(),
            legend.position = "bottom")

ggsave(plot = p3, filename = "./figs/model_results/sfd/summary_att_density_plots.pdf",
       device = "pdf", height = 8.3, width = 11.7, units = "in")

summary_beta_cf <- sim_beta_cf[, quantile(.SD$percent_change, probs = c(0.05, 0.5, 0.95)), by = c("location", "model_name")]

percentiles <- rep(c("q025", "q50", "q975"), dim(args)[1] + dim(args)[1]/4)
summary_beta_cf[, percentile := percentiles]
summary_beta_cf <- dcast(summary_beta_cf, formula = location + model_name ~ percentile, value.var = "V1")

summary_beta_cf[, location := str_to_title(location)]

p4 <- ggplot(summary_beta_cf, aes(x = q50, color = model_name, y = location)) +
                  geom_errorbarh(aes(xmin = q025, xmax = q975), size = 2, height = 0, position = position_dodgev(height = 0.4)) +
                  geom_point(size = 3, position = position_dodgev(height = 0.4)) +
                  geom_vline(xintercept = 0, linetype = 2) +
                  labs(title = "Effect of CFH Policy on Eviction Filings",
                       x = "Percentage change",
                       y = "",
                       color = "Model",
                       caption = "Bars are 95% UI.") +
                  theme_bw() +
                  scale_color_manual(values = c("#4d9221", "#b8e186", "#2166ac", "#92c5de")) +
                  scale_x_continuous(labels = percent, breaks = seq(-0.1, 1, 0.1)) +
                  theme(plot.title = element_text(hjust = 0.5),
                        strip.background = element_blank(),
                        legend.position = "bottom",
                        axis.title = element_text(size = 12),
                        axis.text.y = element_text(size = 12),
                        axis.text.x = element_text(size = 12,
                                                   margin = margin(t = 5, r = 0, b = 10, l = 0)))

ggsave(plot = p4, filename = "./figs/model_results/sfd/att_eviction_filings.pdf",
       device = "pdf", height = 8.3, width = 11.7, units = "in")


############################################################################
# Data analysis:
# What is the difference in eviction rates between block groups w|w/o CFH? #
#############################################################################

df_evict_sum <- copy(df_sfd)

df_evict_sum[, evict_rate_mean := mean(.SD$evict_rate), by = c("location", "cfh_any")]
df_evict_sum[, evict_rate_lower := quantile(.SD$evict_rate, probs = 0.05), by = c("location", "cfh_any")]
df_evict_sum[, evict_rate_upper := quantile(.SD$evict_rate, probs = 0.95), by = c("location", "cfh_any")]

df_evict_sum <- unique(df_evict_sum[, .(location, evict_rate_mean, evict_rate_lower,
                                        evict_rate_upper, cfh_any)])

df_evict_sum[, location := str_to_title(location)]
df_evict_sum[, cfh := ifelse(cfh_any == 1, "t", "f")]

df_evict_sum <- dcast(df_evict_sum, location ~ cfh, value.var = "evict_rate_mean")

plot_colors <- c("cfh_f" = "#a6e6db", "cfh_t" = "#0d2b52")

p5 <- ggplot(df_evict_sum, aes(y = location, color = cfh)) +
        geom_linerange(aes(xmin = f, xmax = t), color = 'black', linewidth = 1, alpha = 0.3) +
        geom_point(size = 4, aes(x = t, color = "cfh_t")) +
        geom_point(size = 4, aes(x = f, color = "cfh_f")) +
        labs(x = "\nEviction Filings per 100 Renting Households",
             y = "",
             title = "Average eviction filings are greater in \nneighborhoods with crime-free housing properties",
             color = "") +
        scale_x_continuous(limits = c(0, 12), breaks = seq(0, 12, 4)) +
        scale_color_manual(values = plot_colors,
                           labels = c("No Crime-Free Housing Properties",
                                      "Neighborhoods with Crime-Free \nHousing Properties")) +
        theme_bw() +
        theme(axis.line.x = element_line(colour = "black", size = 0.65),
              axis.ticks.x = element_line(size = 0.65),
              axis.ticks.length = unit(2, "mm"),
              axis.ticks.y = element_blank(),
              axis.text = element_text(size = 12, family = "sans"),
              axis.title.y = element_text(angle = 0, vjust = 0.5),
              axis.title = element_text(size = 12, family = "sans"),
              legend.position = "bottom",
              legend.justification = "left",
              legend.text = element_text(size = 8, family = 'sans'),
              plot.title = element_text(size = 14, family = 'sans'))

ggsave(plot = p5, filename = "./figs/average_eviction_filings.pdf", 
       device = 'pdf', bg = "transparent",
       height = 8.27/2, width = 11.69/2, unit = 'in')

#######################################
# Maps of evictions and CFH locations #
#######################################

tiger_state  <- st_read("./data/shapefiles/tiger_state.geojson")
tiger_county <- st_read("./data/shapefiles/tiger_county.geojson")
tiger_place  <- st_read("./data/shapefiles/tiger_place.geojson")

plot_map <- function(city){
  
  df_map <- setDT(st_read(paste0("./data/processed/sfd/", city, "_prepped.geojson")))
  df_map[, evict_rate_lab := cut(evict_rate, breaks = c(-Inf, 0, 1, 5, 10, 999), 
                                  labels = c("0", ">0 - 1", ">1 - 5", ">5 - 10", ">10"),
                                  include.lowest = T)]
  
  df_map <- st_as_sf(df_map) %>%
            st_transform(., st_crs(tiger_place))
  
  # Crop map based on tiger file to exact city boundaries:
  tiger_crop <- tiger_place[tiger_place$location == city,]
  
  df_map <- st_crop(df_map, tiger_crop)
  
  # Construct random points for each CFH locations within respective polygons:
  sample_points <- function(geo, s){
    as.data.table(st_sample(geo, size = s, type = "random", exact = T))
  }
  
  df_cfh <- copy(df_map[df_map$cfh_num != 0,])
  
  df_cfh <- mcmapply(sample_points, 
                      geo = df_cfh$geometry, 
                      s = df_cfh$cfh_num,
                      SIMPLIFY = F)
  
  df_cfh <- setDT(rbindlist(df_cfh))
  df_cfh <- st_as_sf(df_cfh) %>%
            st_set_crs(., st_crs(tiger_place))
  
  map_plot <- ggplot(df_map) + 
                geom_sf(aes(fill = evict_rate_lab)) + 
                geom_sf(data = df_cfh, aes(color = cfh_num), size = 2, color = 'black', 
                        shape = 17) +
                theme_void() +
                scale_fill_manual(values = c("#ffffff", "#D4D4D4", "#B4B4B4", "#909090", "#636363")) +
                labs(fill = "Number of Evictions") +
                guides(fill = guide_legend(nrow = 1, byrow = T),
                       color = guide_legend(nrow = 1)) +
                theme(legend.position = "bottom",
                      legend.direction = "vertical",
                      legend.title.align = 0.5,
                      plot.title = element_text(hjust = 0.5, family = 'sans', size = 16),
                      legend.text = element_text(family = 'sans', size = 12),
                      legend.title = element_text(family = 'sans', size = 14))
  
  return(map_plot)
  
}

maps <- lapply(cfh_cities, plot_map)

plots <- ggarrange(plotlist = maps, ncol = 2, nrow = 2, common.legend = T, 
                   legend = "bottom", 
                   labels = str_to_title(cfh_cities),
                   font.label = list(size = 14, family = "sans"))

ggsave("./figs/cfh_maps.pdf", plots, 
       device = "pdf", width = 11.69, height = 8.27/1.4, units = "in")

##
# Location of CFH Cities in Arizona
##

df_treated <- df[year == max(year) & treated == T, .(geoid, treated)]
setnames(df_treated, "geoid", "place_geoid")

tiger_place <- setDT(tiger_place)

df_map <- join(tiger_place, df_treated, by = "place_geoid", type = "left")
df_map <- st_as_sf(df_map[df_map$treated == T,])

tiger_place <- st_as_sf(tiger_place) %>%
               st_set_crs(st_crs(tiger_state))

df_map <- st_transform(df_map, 3488) %>%
            st_centroid(df_map) %>%
            st_buffer(dist = set_units(5, "km")) %>%
            st_transform(st_crs(tiger_state))

df_map$title <- "Cities with a Crime-Free \nHousing Program"

text_city  <- copy(df_map) %>% 
              st_centroid() %>%
              setDT()

# Nudge names of cities
text_city[location == "flagstaff", geometry := geometry + c(-0.3, 0.15)]
text_city[location == "cottonwood", geometry := geometry + c(-0.35, 0.15)]
text_city[location == "peoria", geometry := geometry + c(-0.32, 0.12)]
text_city[location == "scottsdale", geometry := geometry + c(0.37, 0.15)]

text_city[location == "glendale", geometry := geometry + c(-0.35, 0.12)]
text_city[location == "phoenix", geometry := geometry + c(0.35, -0.04)]

text_city[location == "avondale", geometry := geometry + c(-0.39, 0.1)]

text_city[location == "apache junction", geometry := geometry + c(0.65, 0)]
text_city[location == "mesa", geometry := geometry + c(0.48, 0.1)]

text_city[location == "tempe", geometry := geometry + c(-0.33, -0.1)]
text_city[location == "chandler", geometry := geometry + c(-0.32, -0.12)]

text_city[location == "casa grande", geometry := geometry + c(0.5, -0.1)]
text_city[location == "queen creek", geometry := geometry + c(0.45, -0.12)]

text_city[location == "gilbert", geometry := geometry + c(0.6, -0.04)]

text_city[location == "yuma", geometry := geometry + c(0.35, 0.15)]

text_city[location == "marana", geometry := geometry + c(0.37, 0.03)]
text_city[location == "tucson", geometry := geometry + c(-0.35, -0.15)]

text_city[location == "sierra vista", geometry := geometry + c(0.35, 0.15)]

p6  <- ggplot() + 
        geom_sf(data = tiger_state, fill = "#f7f7f7") + 
        geom_sf(data = tiger_county, fill = "transparent", alpha = 0.6, color = "#bdbdbd") + 
        geom_sf(data = df_map, fill = "black", alpha = 0.4) +
        geom_sf_text(data = text_city, aes(label = str_to_title(location), geometry = geometry)) +
        #facet_grid(. ~ title) +
        theme_void() +
        labs(title = "Cities with Crime-Free Housing Programs") +
        guides(fill = guide_legend(nrow = 2, byrow = T)) +
        theme(legend.position = "bottom",
              legend.direction = "vertical",
              legend.title.align = 0.5,
              plot.title = element_text(size = 16, family = "sans", hjust = 0.5))

ggsave("./figs/state_cfh.pdf", plot = p6, 
       device = "pdf", height = 11.69/1.2, width = 8.27/1.2, units = "in")

##
# Time Trends - Number of Renters exposed to the policy
##

df_year <- df[, .(year, location, post_treat, pop_tenants, total_population)]

df_year[, total_pop := sum(.SD$total_population), by = "year"]
df_year <- df_year[post_treat == 1,]

df_year[, cfh_total := sum(.SD$total_population), by = "year"]

df_year <- unique(df_year[, .(cfh_total, year)])

p7 <- ggplot(df_year, aes(x = year, y = cfh_total)) +
                    geom_point(size = 3, color = "black", alpha = 0.7) +
                    geom_line(linewidth = 1.2, color = "black", alpha = 0.7) +
                    scale_y_continuous(labels = label_number(scale = 1e-6),
                                       breaks = seq(0e6, 8e6, 2e6),
                                       limits = c(0e6, 5e6)) +
                    scale_x_continuous(breaks = seq(2009, 2023, 2)) +
                    labs(y = "Population \n(in millions)",
                         x = "Year",
                         title = "21% More Arizonans Live in Crime-Free Housing Cities Since 2009") +
                    theme_bw() +
                    theme(legend.position = "none",
                          axis.ticks = element_line(linewidth = 1),
                          axis.ticks.length = unit(5.6, "points"),
                          axis.text = element_text(family = "sans", size = 12),
                          axis.title = element_text(family = "sans", size = 14),
                          axis.title.y = element_text(angle = 0, vjust = 0.5),
                          title = element_text(size = 16, family = "sans"))

ggsave("./figs/cfh_time_trends_population.pdf", p7, 
       device = "pdf", width = 11.69, height = 8.27/1.4, units = "in")


##################
# Within cities: #
##################

sfd_files <- "./data/processed/sfd/"
files <- list.files(sfd_files)[list.files(sfd_files) %like% ".geojson" &
                               !(list.files(sfd_files) %like% "crime")]
files <- paste0(sfd_files, files)

cfh_cities <- c("avondale", "chandler", "mesa", "scottsdale")

compare_vars <- c("evict_rate", "renter_black", "renter_asian", "renter_white_alone",
                  "renter_native_american", "renter_hispanic_latin", 
                  "median_income_10k", "number_rental_units")

# Use Welch–Satterthwaite approximation to estimate degrees-of-freedom for
# t-dist as a means to derive critical value for estimating CI.

# This is relatively pendantic; critical values are typically in range of
# 2.03 to 1.98 rather than routine critical value from z of 1.96.
welch_satt <- function(sd0, sd1, n0, n1){
  
  d <- ((sd0/n0) + (sd1/n1))^2/(((sd0^2)/((n0^2)*(n0 - 1))) + ((sd1^2)/((n1^2)*(n1 - 1))))
  return(d)
  
}

formatted_results <- function(m, s, ci = F, eff_df = NA, sig = 0.025, t_dist = F, digs = 1){
  
  # If assuming normality, use z-score for critical value (t = F). 
  # If effective degrees-of-freedom are set and t_dist = T, use student t to 
  # determine critical value:
  if (ci == T){
    if (t_dist == T){
      cv <- qt(p = sig, df = eff_df, lower.tail = F)
    }else{
      cv <- qnorm(p = sig, lower.tail = F)
    }
    
    return(paste(round(m, digs), paste0("(", round(m - cv*s, digs), ", ", round(m + cv*s, digs), ")")))
  }else{
    return(paste(round(m, digs), paste0("(", round(s, digs), ")")))
  }
}

calc_sfd_summaries <- function(city){
  
  file <- files[files %like% city]
  df_city <- setDT(st_read(file))
  
  df_city[, evict_rate := (evict_count*100)/number_rental_units]
  df_city <-  melt(df_city, id.vars = c("block_geoid", "cfh_any"), measure.vars = compare_vars)
  
  df_city <- df_city[, `:=`(mean = mean(.SD$value, na.rm = T),
                            sd   = sd(.SD$value, na.rm = T),
                            n    = .N), 
                     by = c("cfh_any", "variable")]
  
  new_names <- c("Eviction Rate per 100 Rental Units", "Renter Pop: % Black",
                 "Renter Pop: % asian", "Renter pop: % white", "Renter pop: % Indigenous",
                 "Renter Pop: % Hispanic/Latin", "Median income", 
                 "Number of Rental Units")
  
  df_city[, variable := factor(mapvalues(variable, unique(df_city$variable), new_names), levels = new_names)]
  df_city <- unique(df_city[, .(cfh_any, variable, mean, sd, n)])
  
  # Convert proportions to percentages:
  df_city[variable %like% "%", `:=`(mean = mean*100, sd = sd*100)]
  
  # Convert median income into raw amounts:
  df_city[variable %like% "income", `:=`(mean = mean*1e4, sd = sd*1e4)]
  
  control_n <- unique(df_city[cfh_any == 0]$n)
  treat_n   <- unique(df_city[cfh_any == 1]$n)
  
  df_city <- dcast(df_city, variable ~ cfh_any, value.var = c("mean", "sd", "n"))
  
  df_city[, mean_diff := mean_1 - mean_0]
  
  df_city[, sd_diff := sqrt(((sd_0^2)/n_0) + ((sd_1^2)/n_1))]
  df_city[, df_diff := welch_satt(sd_0, sd_1, n_0, n_1)]
  
  df_city[, Control := formatted_results(mean_0, sd_0)]
  df_city[, Treated := formatted_results(mean_1, sd_1)]
  df_city[, Difference := formatted_results(mean_diff, sd_diff, ci = T, eff_df = df_diff)]
  
  df_city <- df_city[, c("variable", "Control", "Treated", "Difference"), with = F]
  
  df_city <- rbindlist(list(df_city, list("N", control_n, treat_n, "")), use.names = F)
  
  df_city[, City := city]
  
  return(df_city)
  
}

sfd_summary <- rbindlist(lapply(cfh_cities, calc_sfd_summaries))
sfd_summary <- sfd_summary[, .(City, variable, Control, Treated, Difference)]

sfd_summary[, City := mapvalues(City, unique(sfd_summary$City),
                                c("Avondale", "Chandler", "Mesa", "Scottsdale"))]

write.csv(sfd_summary, "./figs/city_summaries_2023.csv", row.names = F)


#######################################
# Demographic differences by location #
#######################################

df <- fread("./data/processed/df_crime_analysis.csv")

df_z <- copy(df[year == 2023, .(location, treated, pop_black, pop_asian, pop_latin_hispanic, pop_white, total_population)])

df_z <- melt(df_z, 
             id.vars = c("location", "treated", "total_population"), 
             value.name = "prop",
             variable.names = "var")

df_z[, treated := ifelse(treated == T, 1, 0)]
df_z[, variable := factor(variable, levels = rev(c("pop_black", "pop_asian", "pop_latin_hispanic", 
                                                   "pop_white")))]

# Then calculate median and lower/upper notches:
df_z[, mid := mean(.SD$prop), by = c("variable", "treated")]
df_z[, lower := t.test(.SD$prop)$conf.int[1], by = c("variable", "treated")]
df_z[, upper := t.test(.SD$prop)$conf.int[2], by = c("variable", "treated")]

plot_colors = c("#ba9c7d", "#807dba")

# Remap variables:
df_z[, treated := ifelse(treated == 1, "Yes", "No")]

new_names <- c("Black", "Asian", 
               "Hispanic", "White")

df_z[, variable := factor(mapvalues(variable, unique(df_z$variable), new_names), levels = new_names)]

pop_prop <- ggplot(df_z, aes(x = prop, y = treated, color = treated)) +
  facet_wrap(~variable,  nrow = 4, strip.position = 'right') +
  geom_point(alpha = 0.3, position = position_jitterdodge()) +
  geom_crossbar(aes(xmin = lower, x = mid, xmax = upper)) +
  theme_bw() +
  labs(title = "Cities with crime-free housing programs have larger Black and Asian populations",
       x = "Proportion of total population",
       y = "City with CFHP") +
  scale_color_manual(values = plot_colors) +
  scale_x_continuous(breaks = seq(0, 1, 0.2), labels = percent) +
  theme(strip.text.y.right = element_text(size = 14, family = 'serif', angle = 0),
        strip.background = element_blank(),
        axis.ticks = element_line(linewidth = 1),
        axis.ticks.length = unit(5.6, "points"),
        axis.text = element_text(family = "serif", size = 12),
        axis.title = element_text(family = "serif", size = 18),
        axis.title.y = element_text(angle = 0, vjust = 1),
        panel.spacing = unit(0, "lines"),
        legend.position = "none")

ggsave(pop_prop, filename = "./figs/demographic_differences_2023.pdf", 
       device = "pdf", bg = "transparent", width = 11.7, height = 8.3, units = "in")
