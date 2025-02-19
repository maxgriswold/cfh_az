library(conley)
library(data.table)
source("./code/construct_sfd.R")

city_names <- c()

run_sfd_models <- function(city_name, dep, ind, model_name){
  
  sfd_save <- paste0("./cfho_analysis/sfd/", city_name, "_prepped_sfd.geojson")
  df_city <- st_read(sfd_save)
  
  # Check sensitivity of model to targeted &  reduced controls - do we arrive
  # at similar findings, if we limit to only immediate neighbors of treated sites 
  # for constructing comparisons?
  
  if (analysis_type == "sensitivity_test_neighbors"){
    
    df_city['index'] <- 1:dim(df_city)[1]
    treat_sites <- df_city[df_city$cfho_any == 1,]$index
    
    neigh <- st_intersects(df_city, df_city)
    neigh <- neigh[treat_sites]
    neigh <- unique(c(unlist(neigh), treat_sites))
    
    df_city <- df_city[df_city$index %in% neigh,]
    
  }
  
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
  
  spec <- paste0(paste0("sfd.", dep), " ~ ", paste0("sfd.", ind, collapse = " + "))
  
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