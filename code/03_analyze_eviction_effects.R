library(conleyreg)
library(data.table)
library(sf)

source("./code/construct_sfd.R")

model_dir <- "./data/models"

run_sfd_models <- function(city_name, dep, ind, model_name, model_type = "sfd"){
  
  sfd_save <- paste0("./data/processed/sfd/", city_name, "_prepped.geojson")
  df_city <- st_read(sfd_save)
  
  # Hold onto minimum variables:
  df_city <- df_city[, c(dep, ind, "geometry", "block_geoid")]
  
  sfd_city <- construct_sfd(spatial_df = df_city,
                            dependent_var = dep,
                            independent_vars = ind,
                            rotation = seq(0, 330, 30),
                            plot = F,
                            tol = 0)
  
  # Set up summary results table which holds overall estimates across models
  # and individual model results for a given angle:
  full_results <- list()
  full_results[['sfd_data']] <- sfd_city
  
  # For each angle, use calculated spatial first-differences within a OLS
  # regression using conley standard errors (i.e. account for clustered spatial
  # correlations)
  
  # Calculate centroid for each block, then add as lat/lon to the original dataset:
  df_conley <- st_centroid(df_city)
  
  # Note: We need to index rows in the original dataset to merge w/ the SFD
  # calculations:
  setDT(df_conley)
  
  df_conley[, lon := st_coordinates(geometry)[, 1]]
  df_conley[, lat := st_coordinates(geometry)[, 2]]
  
  df_conley[, rowno := .I]
  
  df_conley <- df_conley[, .(block_geoid, lon, lat, rowno)]
  
  # Add on location
  sfd_city <- join(sfd_city, df_conley, by = "rowno", type = "left")
  
  if (model_type %in% c("sfd", "sdd")){
    spec <- paste0(paste0(paste0(model_type, "."), dep), " ~ ", paste0(paste0(model_type, "."), ind, collapse = " + "))
  }else{
    stop("Error: model_type must be either sfd or sdd")
  }
  
  for (a in unique(sfd_city$angle)){
    
    df_a <- sfd_city[sfd_city$angle == a,]
    
    # Get dataframe of coefficient table
    lm_con <- conleyreg(spec, data = df_a, lat = 'lat', lon = 'lon', dist_cutoff = 200, 
                        gof = T, vcov = F)
    
    coef_res <- as.data.table(lm_con$coef[,])
    coef_res[, var := c("intercept", ind)]
    coef_res[, angle := a]
    coef_res[, city := city_name]
    
    sum_res <- data.table("nobs" = lm_con$nobs,
                          "rsquared" = lm_con$r.squared,
                          "rsquared_adj" = lm_con$adj.r.squared,
                          "model_name" = model_name)
    
    # Get dataframe of covariance-variance matrix
    lm_con <- conleyreg(spec, data = df_a, lat = 'lat', lon = 'lon', dist_cutoff = 200, 
                        gof = T, vcov = T)
    
    vc_res <- as.data.table(lm_con$vcov)
    vc_res[, var := c("intercept", ind)]
    vc_res[, angle := a]
    vc_res[, city := city_name]
    
    res_list = list("coef" = coef_res, "vcov" = vc_res, "sum" = sum_res)
    
    results_name <- paste0("angle=", a)
    full_results[[results_name]] <- res_list
    
  }
  
  return(full_results)
  
}


# Model settings

cfh_cities    <- c("avondale", "chandler", "mesa")
dep_var       <- "evict_rate"
adjusted_spec <- c("median_income_10k", "renter_white_alone", "renter_black", "renter_asian",
                   "renter_hispanic_latin", "percent_poverty_150")

#
# Binary, Unadjusted
#

sfd_results_any_unadj <- lapply(cfh_cities, run_sfd_models, 
                                dep = dep_var, 
                                ind = c("cfh_any"),
                                model_name = "cfh_any_unadj")

names(sfd_results_any_unadj) <- cfh_cities
saveRDS(sfd_results_any_unadj, sprintf("%s/cfh_any_unadj.rds", model_dir))
        
#
# Binary, adjusted
#

sfd_results_any_adj   <- lapply(cfh_cities, run_sfd_models, 
                                dep = dep_var, 
                                ind = c("cfh_any", adjusted_spec),
                                model_name = "cfh_any_adj")

names(sfd_results_any_adj) <- cfh_cities
saveRDS(sfd_results_any_adj, sprintf("%s/cfh_any_adj.rds", model_dir))
        
#
# Continuous, unadjusted
#
      
sfd_results_count_unadj <- lapply(cfh_cities, run_sfd_models, 
                                  dep = dep_var,
                                  ind = c("cfh_num"),
                                  model_name = "cfh_count_unadj")

names(sfd_results_count_unadj) <- cfh_cities
saveRDS(sfd_results_count_unadj, sprintf("%s/cfh_count_unadj.rds", model_dir))

#
# Continuous, adjusted
#

sfd_results_count_adj   <- lapply(cfh_cities, run_sfd_models, 
                                  dep = dep_var,
                                  ind = c("cfh_num", adjusted_spec),
                                  model_name = "cfh_count_adj")

names(sfd_results_count_adj) <- cfh_cities
saveRDS(sfd_results_count_adj, sprintf("%s/cfh_count_adj.rds", model_dir))
        
