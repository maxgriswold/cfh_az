rm(list = ls())

library(data.table)
library(ggplot2)
library(ggstance)
library(ggridges)
library(MASS)
library(stringr)
library(scales)

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

df_sfd <- setDT(ldply(sfd_datasets, combine_sfd_df))
df_sfd[, geometry := NULL]

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

args <- expand.grid(loc = unique(sim_results$location),
                    modelname = unique(sim_results$model_name))

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

ggsave(plot = p2, filename = "./figs/model_results/sfd/summary_att_density_plots.pdf",
       device = "pdf", height = 8.3, width = 11.7, units = "in")

summary_beta_cf <- sim_beta_cf[, quantile(.SD$percent_change, probs = c(0.05, 0.5, 0.95)), by = c("location", "model_name")]

percentiles <- rep(c("q025", "q50", "q975"), dim(args)[1] + dim(args)[1]/3)
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


##############################################################################################
# Data analysis:
# What is the difference in eviction rates between block groups with and without the policy? #
##############################################################################################



