# Use augsynth to determine if crime-free housing led to reductions in crime
# Max Griswold
# 2/10/25

rm(list = ls())

library(augsynth)
library(data.table)
library(ggplot2)
library(ggrepel)
library(kableExtra)

library(did)
library(plyr)
library(gridExtra)

library(modelsummary)
library(stargazer)

# Code options:
remove_outliers   <- T
run_placebo       <- F
run_lambda_search <- F
estimate_did      <- T

df <- fread("./data/processed/df_crime_analysis.csv")

# Removealways-treated locations, and single pre-treatment locations
df <- df[implementation_date > 2010|is.na(implementation_date),]

# Sensitivity test: Does removing outliers w/ implausible UCR data change
# estimated results?

if (remove_outliers){
  df <- df[outlier == F,]
}

# Set up model options
outcome_vars <- c("index_total", "assault_total", "burglary_total")
covs         <- c("pop_white", "pop_black", "pop_asian", "pop_native_american",
                  "pop_latin_hispanic", "median_income_adj", "total_population")

# There is data missing for a single row in the median_income variable.
# within the city of Quartzsite. While I could use multiple imputation
# (e.g., MICE) to impute this value, median income is relatively stable
# during this time. So, I am putting a linearly interpolated value for this
# missing row to ensure a complete time series across all units. I did not
# find any meaningful difference in results when I commented this row of code out.

df[location == "Quartzsite" & year == 2021, median_income_adj := 16138.65]

# Wrapper code for multisynth
run_multi <- function(out, dd, covars = NULL, l = 0){
  
  if (is.null(covars)){
    spec <- formula(paste0(out, " ~ post_treat"))
  }else{
    spec <- formula(paste0(out, " ~ post_treat |", paste0(covars, collapse = " + ")))
  }
  
  mod <- multisynth(form = spec, unit = location, time = year,
                    lambda = l, data = dd)
  
  res <- summary(mod)
  
  return(res)
  
}

# Models w/o covariates
args <- data.frame("out" = outcome_vars)
mods_1 <- mlply(args, run_multi, dd = df)
names(mods_1) <- outcome_vars

# Models w/ covariates
args <- data.frame("out" = outcome_vars)
mods_2 <- mlply(args, run_multi, covars = covs, dd = df)
names(mods_2) <- outcome_vars

# Results appear to be relatively stable over initial set of results,
# regardless of covariate inclusion. 

# Based on these initial findings,let's try out a series of sensitivity tests. 
# First, let's implement a placebo tests to determine
# if the estimated effect changes substantially if a control groups had been "treated".

# Then, see if model fit can be improved through reducing imbalance 
# in pretreatment byfitting the model across a range of lambda values (more/less regularization).

# Ben-Michael, E., Feller, A., & Rothstein, J. (2021). The Augmented Synthetic Control Method. 
# Journal of the American Statistical Association, 116(536), 1789–1803. 
# https://doi.org/10.1080/01621459.2021.1929245

# Pickett, R. E., Hill, J., & Cowan, S. K. (2022). The Myths of Synthetic Control: Recommendations for Practice.
# https://as.nyu.edu/content/dam/nyu-as/cashtransferlab/documents/RPJHSC_SyntheticControl_041122.pdf

################# 
# Placebo tests #
#################

# For each control unit, set each unit to treated, then
# re-estimate the model. Hold onto ATT and se, then plot the range of results.

# Alternative tests to consider:
# Placebo treatment_years: Set treatment to earlier for treated units
# Placebo control w/ year: randomize implementation year and unit

if (run_placebo){
  
  placebo_runs <- function(out, dd, covars = NULL){
    
    control_units <- unique(dd$location[dd$treated == F])
    
    placebo_res <- list()
    
    for (cu in control_units){
      
      dd_c <- copy(dd)
      
      # For this test, I'm setting the earliest year of implementation
      # within treated units as the implementation year. I'm doing this
      # as a "conservative" test since this leads to less possible
      # balancing of the placebo.
      
      dd_c[location == cu, `:=`(treated = T, implementation_date = 2014,
                                post_treat = c(rep(0, 5), rep(1, 10)))]
      
      # Adding small lambda value to push past non-convexity
      # Troublesome model: burglary_total w/ Wickenburg as placebo unit 
      model <- run_multi(out = out, dd_c, covars = covars, l = 0.1)
      effs <- model$att
      
      # Get average post-treatment effect w/ 95% UI:
      res <- data.table(outcome = out,
                        location = cu,
                        covariates = ifelse(is.null(covars), F, T),
                        mean_effect = effs[effs$Level=="Average" & is.na(effs$Time),]$Estimate,
                        lower = effs[effs$Level=="Average" & is.na(effs$Time),]$lower_bound,
                        upper = effs[effs$Level=="Average" & is.na(effs$Time),]$upper_bound)
      
      placebo_res[[cu]] <- res
      
    }
    
    # Add on baseline model results:
    model <- run_multi(out = out, dd, covars = covars)
    effs <- model$att
    
    res <- data.table(outcome = out,
                      location = "baseline",
                      covariates = ifelse(is.null(covars), F, T),
                      mean_effect = effs[effs$Level=="Average" & is.na(effs$Time),]$Estimate,
                      lower = effs[effs$Level=="Average" & is.na(effs$Time),]$lower_bound,
                      upper = effs[effs$Level=="Average" & is.na(effs$Time),]$upper_bound)
    
    placebo_res[["baseline"]] <- res
    
    placebo_res <- rbindlist(placebo_res)
    
    return(placebo_res)
    
  }
  
  plot_placebo <- function(dd){
    
    baseline <- dd[location == "baseline"]
    outcome  <- unique(dd$outcome)
    
    cov_name <- ifelse(unique(dd$covars) == T, ", with covariates", "")
    plot_title <- paste0("Placebo Test: ", outcome, cov_name)
    
    dd <- dd[location != "baseline",]
    
    plot <- ggplot(dd, aes(y = location, x = mean_effect, xmin = lower, xmax = upper)) +
              geom_linerange(alpha = 0.9) +
              geom_point(size = 2) +
              geom_vline(linetype = 2, xintercept = baseline$mean_effect, color = 'red') +
              geom_vline(linetype = 2, xintercept = baseline$lower, color = 'red') +
              geom_vline(linetype = 2, xintercept = baseline$upper, color = 'red') +
              labs(x = "Estimated Average ATT",
                   y = "Placebo Unit",
                   title = plot_title)
              theme_bw() 
    
    return(plot)
    
  }
  
  # Models w/o covariates
  args <- data.frame("out" = outcome_vars)
  placebo_1 <- mlply(args, placebo_runs, dd = df)
  
  # Models w/ covariates
  args <- data.frame("out" = outcome_vars)
  placebo_2 <- mlply(args, placebo_runs, covars = covs, dd = df)
  
  # Plot results for each model
  all_placebos <- c(placebo_1, placebo_2)
  plots <- lapply(all_placebos, plot_placebo)
  
  fig_title <- ifelse(remove_outliers, "placebo_tests_outliers_removed", "placebo_tests_full_data")
  
  ggsave(sprintf("./figs/model_results/multisynth/%s.pdf", fig_title), 
         plot = marrangeGrob(plots, nrow = 1, ncol = 1, top = NULL), 
         height = 11.69, width = 8.27)
  
}

# Placebo tests w/ and w/o outliers are really encouraging - none of the models
# lead to radically different estimates across tests. Accordingly, moving on
# to improving balance:

#################
# Lambda search #
#################

if (run_lambda_search){
  
  lambda_search <- function(out, dd, ls, covars = NULL){
    
    res <- list()
    
    if (is.null(covars)){
      spec <- formula(paste0(out, " ~ post_treat"))
    }else{
      spec <- formula(paste0(out, " ~ post_treat |", paste0(covars, collapse = " + ")))
    }
    
    for (l in ls){
      
      mod <- multisynth(spec, unit = location, time = year,
                        data = dd, lambda = l)
      
      effs <- summary(mod)$att
      imbalance <- mod$global_l2
      
      item_name <- paste(l)
      res[[item_name]] <- data.table("outcome" = out, "covariates" = ifelse(is.null(covars), F, T),
                                     "lambda" = l, "L2 Imbalance" = imbalance, 
                                     "mean_effect" = effs[effs$Level=="Average" & is.na(effs$Time),]$Estimate,
                                     "lower" = effs[effs$Level=="Average" & is.na(effs$Time),]$lower_bound,
                                     "upper" = effs[effs$Level=="Average" & is.na(effs$Time),]$upper_bound)
    }
    
    res <- rbindlist(res)
    return(res)
    
  }
  
  lambda_plot <- function(dd){
    
    outcome  <- unique(dd$outcome)
    
    cov_name <- ifelse(unique(dd$covariates) == T, ", with covariates", "")
    plot_title <- paste0("Lambda search: ", outcome, cov_name)
    
    plot <- ggplot(dd, aes(x = lambda, y = `L2 Imbalance`)) +
              geom_point() +
              labs(x = "Lambda", y = "L2 Imabalance",
                   title = plot_title) +
              theme_bw()
              
    return(plot)
              
  }
  
  lambda_effects <- function(dd){
    
    outcome  <- unique(dd$outcome)
    
    cov_name <- ifelse(unique(dd$covariates) == T, ", with covariates", "")
    plot_title <- paste0("Lambda search: ", outcome, cov_name)
    
    dd[, lambda := as.factor(lambda)]
    
    plot <- ggplot(dd, aes(y = lambda, x = mean_effect, xmin = lower, xmax = upper)) +
              geom_linerange(alpha = 0.9) +
              geom_point(size = 2) +
              labs(x = "Estimated Average ATT",
                   y = "Lambda",
                   title = plot_title)
              theme_bw() 
    
    return(plot)
    
  }
  
  lambda_vals <- c(0, 0.1, 0.2, 0.5, 1, 2, seq(5, 45, 5), seq(50, 100, 10))
  
  # Models w/o covariates
  args <- data.frame("out" = outcome_vars)
  lambda_1 <- mlply(args, lambda_search, dd = df, ls = lambda_vals)
  
  # Models w/ covariates
  args <- data.frame("out" = outcome_vars)
  lambda_2 <- mlply(args, lambda_search, covars = covs, dd = df, ls = lambda_vals)
  
  # Plot results for each model
  all_lambdas <- c(lambda_1, lambda_2)
  plots <- lapply(all_lambdas, lambda_plot)
  
  fig_title <- ifelse(remove_outliers, "lambda_search_outliers_removed", "lambda_search_full_data")
  
  ggsave(sprintf("./figs/model_results/multisynth/%s.pdf", fig_title), 
         plot = marrangeGrob(plots, nrow = 1, ncol = 1, top = NULL), 
         width = 11.69, height = 8.27)
  
  plots <- lapply(all_lambdas, lambda_effects)
  fig_title <- ifelse(remove_outliers, "lambda_v_effects_outliers_removed", "lambda_v_effects_full_data")
  
  ggsave(sprintf("./figs/model_results/multisynth/%s.pdf", fig_title), 
         plot = marrangeGrob(plots, nrow = 1, ncol = 1, top = NULL), 
         height = 11.69, width = 8.27)
  
}

# Based on plots, unsurprisingly, lambda should be zero for models without covariates.
# Estimated effects do not appear to differ substantially based on lambda value

# with covariates, it looks like 

# outliers removed:
# index total -> 75
# assault_total -> 80
# burglary_total -> 90

# full data:
# index total -> 100
# assault_total -> 100
# burglary_total -> 1

################
# Final models #
################

if (remove_outliers){
  lambda_vals <- list("index_total" = 75, "assault_total" = 80, "burglary_total" = 100)
}else{
  lambda_vals <- list("index_total" = 100, "assault_total" = 100, "burglary_total" = 1)
}

df_att <- list()
model_tables <- list()

for (o in outcome_vars){
  
  mod <- run_multi(o, df, covars = covs, l = lambda_vals[[o]])
  res <- mod$att

  model_tables[[o]] <- mod
  
  # Get average of crime rates in 2023 to calculate percentage effects
  crime_obs <- df[year == 2023, mean(get(o)), ]
  
  df_res <- data.table("outcome" = o,
                      "avg_obs" = crime_obs,
                      "mean_effect" = res[res$Level=="Average" & is.na(res$Time),]$Estimate,
                      "lower" = res[res$Level=="Average" & is.na(res$Time),]$lower_bound,
                      "upper" = res[res$Level=="Average" & is.na(res$Time),]$upper_bound)
  
  df_res[, `:=`(percent_mean = (avg_obs + mean_effect)/avg_obs - 1,
                percent_lower = (avg_obs + lower)/avg_obs - 1,
                percent_upper = (avg_obs + upper)/avg_obs - 1)]
  
  df_att[[o]] <- df_res
  
}

df_att <- rbindlist(df_att)

df_att[, outcome_nice := c("Total Crime", "Assaults", "Burglaries")]

# Plot balance and estimated effects:
plot_bal <- function(o){
  
  mod <- model_tables[[o]]
  
  # Hold onto lags for minimum across the treated groups:
  ests <- mod$att
  setDT(ests)
  ests <- ests[Time >= -5,]
  ests[, label := ifelse(Time == 3, Level, NA),]
  
  ests[, is_avg := ifelse(Level == "Average", "A", "B")]
  
  o_name <- ifelse(o == "index_total", "total crime", ifelse(o == "assault_total", "total assault", "total burglary"))
  
  # Modified plot from the augsynth package
  bal <-  ggplot(ests, aes(x = Time, y = Estimate, group = Level,
                   color = is_avg, alpha = is_avg)) +
        geom_point(size = 1) +
        geom_line(size = 1) +
        geom_label_repel(aes(label = label), nudge_x = 1.3, na.rm = T) +
        geom_hline(yintercept = 0, lty = 2) +
        geom_vline(xintercept = 0, lty = 2) +
        geom_ribbon(data = ests[Level == "Average",], aes(ymin = lower_bound, ymax = upper_bound), alpha = 0.2, fill = "black", color = NA) +
        scale_alpha_manual(values = c(1, 0.5)) +
        scale_color_manual(values=c("#333333", "#818181")) +
        scale_x_continuous(limits = c(-5, 5),
                           breaks = seq(-5, 3, 2),
                           expand = c(0, 2)) +
        guides(alpha = F, color = F) +
        theme_bw() +
        labs(x = "Time relative to treatment",
             y = "Estimate",
             title = sprintf("Estimated effect of CFH on %s rates per 1k people", o_name)) +
        theme(axis.title.y = element_text(angle = 0, vjust = 0.5))
  
  return(bal)
  
}

plots <- lapply(outcome_vars, plot_bal)

ggsave("./figs/model_results/multisynth/balance_plots.pdf",
       plot = marrangeGrob(plots, nrow = 1, ncol = 1, top = NULL),
       width = 11.69,
       height = 8.27)

# Save table results
combine_ests <- function(o){
  
  res <- model_tables[[o]]$att
  res <- res[Level == "Average",]
  res[, outcome := o]
  
  return(res)
  
}

# Prep result table for methods:

res_table <- rbindlist(lapply(outcome_vars, combine_ests))
res_table <- res_table[Time >= 0|is.na(Time),]

# For overall, set Event time == 4 so it appears last in the table
# Recode to average effect later:
res_table[is.na(Time), Time := 4]

res_table <- res_table[, .(outcome, Time, Estimate, Std.Error)]
setnames(res_table, names(res_table), c("outcome", 'event_time', "est", "se"))

res_table[, pval := 2*(1 - pnorm(abs(est/se)))]
res_table[, stars := ifelse(pval < 0.01, "$^{***}$", ifelse(pval < 0.05, "$^{**}$", ifelse(pval < 0.1, "$^{*}$", "")))]
res_table[, est := paste0(round(est, 2), stars, "\n(", round(se, 2), ")")]

res_table[, est := linebreak(est, align = "c")]

res_table <- dcast(res_table, event_time~outcome, value.var = "est")

res_table[, event_time := as.character(event_time)]
res_table[event_time == "4", event_time := "Average"]

res_table[event_time == "0", event_time := "Time Period = 0"]
res_table <- res_table[, .(event_time, index_total, assault_total, burglary_total)]
  
res <- kable(res_table, format = 'latex', booktabs = T, escape = F,
             col.names = c("", "Total Crimes", "Total Assaults", "Total Burglaries"), align = c("rccc"),
             linesep = "") %>%
        row_spec(1, extra_latex_after = "\\addlinespace") %>%
        row_spec(2, extra_latex_after = "\\addlinespace") %>%
        row_spec(3, extra_latex_after = "\\addlinespace") %>%
        row_spec(4, extra_latex_after = "\\addlinespace")

sink("./figs/model_results/multisynth/result_table.txt")
cat(res)
sink()

p <- ggplot(df_att, aes(y = outcome_nice, x = percent_mean, xmin = percent_lower, xmax = percent_upper)) +
        geom_point(size = 3) +
        geom_linerange(linewidth = 1, alpha = 0.3) +
        geom_vline(xintercept = 0, linetype = 2) +
        scale_x_continuous(labels = scales::percent, limits = c(-1, 1.2)) +
        labs(x = "Percent change in crime \nrate per 1,000 people in 2023",
             y = "",
             title = "Crime-Free Housing has no meaningful effect on crime in Arizona") +
        theme_bw() +
        theme(axis.text = element_text(family = 'sans', size = 10),
              plot.title = element_text(family = 'sans', size = 14),
              axis.title = element_text(family = 'sans', size = 12))

ggsave(plot = p, filename = "./figs/estimated_cfh_effects_crime.pdf", height = 8.27/2, width = 8.27)

# Additional check: Does using DID (Callaway and Sant'anna) lead to substantially different estimates? 

if (estimate_did){
  
  df[is.na(implementation_date), implementation_date := 0]
  df[, location := as.numeric(as.factor(location))]
  
  csa_total_unadj <- att_gt(yname = "index_total", 
                      tname = "year", 
                      idname = "location", 
                      gname = "implementation_date", 
                      control_group = "notyettreated", 
                      data = df, 
                      xformla = NULL,
                      est_method = "dr", 
                      bstrap = T,
                      biters = 1000,
                      base_period	= "universal",
                      allow_unbalanced_panel = F)
  
  att_total_unadj <- aggte(csa_total_unadj, type = 'dynamic', balance_e = 3, 
                           na.rm = T, min_e = -5)
  
  csa_assault_unadj <- att_gt(yname = "assault_total", 
                        tname = "year", 
                        idname = "location", 
                        gname = "implementation_date", 
                        control_group = "notyettreated", 
                        data = df, 
                        xformla = NULL,
                        est_method = "dr", 
                        bstrap = T,
                        biters = 1000,
                        base_period	= "universal",
                        allow_unbalanced_panel = F)
  
  att_assault_unadj <- aggte(csa_assault_unadj, type = 'dynamic', balance_e = 3, 
                             na.rm = T, min_e = -5)
  
  csa_burglarly_unadj <- att_gt(yname = "burglary_total", 
                          tname = "year", 
                          idname = "location", 
                          gname = "implementation_date", 
                          control_group = "notyettreated", 
                          data = df, 
                          xformla = NULL,
                          est_method = "dr", 
                          bstrap = T,
                          biters = 1000,
                          base_period	= "universal",
                          allow_unbalanced_panel = F)
  
  att_burglarly_unadj <- aggte(csa_burglarly_unadj, type = 'dynamic', balance_e = 3, 
                               na.rm = T, min_e = -5)
  
  # Add on city-level covariates.
  cov_form <- formula(paste("~ ", paste(covs, collapse = " + ")))
  
  csa_total <- att_gt(yname = "index_total", 
                  tname = "year", 
                  idname = "location", 
                  gname = "implementation_date", 
                  control_group = "notyettreated", 
                  data = df, 
                  xformla = cov_form,
                  est_method = "dr", 
                  bstrap = T,
                  biters = 1000,
                  base_period	= "universal",
                  allow_unbalanced_panel = F,
                  clustervars="location")
  
  att_total <- aggte(csa_total, type = 'dynamic', balance_e = 3, 
                     na.rm = T, min_e = -5)
  
  csa_assault <- att_gt(yname = "assault_total", 
                  tname = "year", 
                  idname = "location", 
                  gname = "implementation_date", 
                  control_group = "notyettreated", 
                  data = df, 
                  xformla = cov_form,
                  est_method = "dr", 
                  bstrap = T,
                  biters = 1000,
                  base_period	= "universal",
                  allow_unbalanced_panel = F)
  
  att_assault <- aggte(csa_assault, type = 'dynamic', balance_e = 3, 
                       na.rm = T, min_e = -5)
  
  csa_burglarly <- att_gt(yname = "burglary_total", 
                        tname = "year", 
                        idname = "location", 
                        gname = "implementation_date", 
                        control_group = "notyettreated", 
                        data = df, 
                        xformla = cov_form,
                        est_method = "dr", 
                        bstrap = T,
                        biters = 1000,
                        base_period	= "universal",
                        allow_unbalanced_panel = F)
  
  att_burglarly <- aggte(csa_burglarly, type = 'dynamic', balance_e = 3, 
                         na.rm = T, min_e = -5)
  
  collect_mods <- list(att_total_unadj, att_total, 
                       att_assault_unadj, att_assault,
                       att_burglarly_unadj, att_burglarly)
  
  model_types <- c("unadjusted", "adjusted") 
  outcomes <- c("total", "assault", "burglarly")
  
  names(collect_mods) <- as.vector(outer(model_types, outcomes, FUN = paste, sep = "_"))
  
  event_table <- function(mods, titl){
    
    res <- list()
    
    for (i in 1:2){
      
      mod <- mods[[i]]
      
      pval <- 2 * (1 - pnorm(abs(mod$att.egt/mod$se.egt)))
      stars <- ifelse(pval < 0.01, "$^{***}$", ifelse(pval < 0.05, "$^{**}$", ifelse(pval < 0.1, "$^{*}$", "")))
      est <- paste0(round(mod$att.egt, 3), stars, "\n(", round(mod$se.egt, 3), ")")
      
      # Add on overall ATT as well:
      pval <- 2 * (1 - pnorm(abs(mod$overall.att/mod$overall.se)))
      stars <- ifelse(pval < 0.01, "$^{***}$", ifelse(pval < 0.05, "$^{**}$", ifelse(pval < 0.1, "$^{*}$", "")))
      est_avg<- paste0(round(mod$overall.att, 3), stars, "\n(", round(mod$overall.se, 3), ")")
      
      res[[i]] <- c(est, est_avg)
      
    }
    
    dd <- data.table("Time Period" = c(mod$egt, 4),
                     Unadjusted = res[[1]],
                     Adjusted = res[[2]])
    
    dd <- dd[`Time Period` >= 0 & `Time Period` <= 4,]
    dd[, `Time Period` := as.character(`Time Period`)]
    dd[`Time Period`== "0", `Time Period` := "Time Period = 0"]
    
    # Note that period "4" is actually the average effect extracted earlier:
    dd[`Time Period`== "4", `Time Period` := "Average"]
    
    
    clusts <- mods[[1]]$DIDparams$n
    obs    <- dim(mods[[1]]$DIDparams$data)[[1]]
    
    dd_add <- data.table("Time Period" = c("Observations", "Clusters", "Adjusted"),
                         "Unadjusted" = c(obs, clusts, "No"),
                         "Adjusted" = c(obs, clusts, "Yes"))
    
    dd <- rbind(dd, dd_add)
    
    dd$Adjusted <- linebreak(dd$Adjusted, align = "c")
    dd$Unadjusted <- linebreak(dd$Unadjusted, align = "c")
    
    res <- kable(dd, format = 'latex', booktabs = T, escape = F,
                 col.names = c("", "(1)", "(2)"), align = c("rcc"),
                 linesep = "") %>%
      row_spec(1, extra_latex_after = "\\addlinespace") %>%
      row_spec(2, extra_latex_after = "\\addlinespace") %>%
      row_spec(3, extra_latex_after = "\\addlinespace") %>%
      row_spec(4, extra_latex_after = "\\addlinespace") %>%
      row_spec(5, extra_latex_after = "\\addlinespace\\midrule")
    
    return(res)
    
  }
  
  table_total    <- event_table(list(att_total_unadj, att_total))
  table_assault  <- event_table(list(att_assault_unadj, att_assault))
  table_burglary <- event_table(list(att_burglarly_unadj, att_burglarly))
  
  sink("./figs/model_results/csa/result_tables.txt")
  cat("Total Crime: 
    ")
  cat(table_total)
  cat("

Total Assaults: 
    ")
  cat(table_assault)
  cat("

Total Burglaries: 
    ")
  cat(table_burglary)
  sink()
  
  # Combine event-time estimates into single dataframe to make
  # plot for checking parallel trends. 
  
  extract_ests <- function(out, mod_type){
    
    mod <- paste0(mod_type, "_", out)
    mod <- collect_mods[[mod]]
    
    dd <- data.table(event_time = mod$egt,
                     est = mod$att.egt,
                     se = mod$se.egt,
                     outcome = out,
                     model = mod_type)
    
    dd[, ci_lower := est - 1.96*se]
    dd[, ci_upper := est + 1.96*se]
    
    return(dd)
    
  }
  
  args <- expand.grid(out = outcomes, mod_type = model_types)
  
  dd_history <- mapply(extract_ests, out = args$out, mod_type = args$mod_type, SIMPLIFY = F)
  dd_history <- rbindlist(dd_history)
  
  # Remove rows for the reference period and outside the range of
  # CSA bounds
  dd_history <- dd_history[event_time != -1 & event_time <= 5 & event_time >= -5,]
  dd_history[, outcome := mapvalues(outcome, from = c("total", "assault", "burglarly"), 
                                    to = c("Total crime rate", "Total assault rate","Total burglarly rate"))]
  
  dd_history[, model := mapvalues(model, from = c("unadjusted", "adjusted"), to = c("CSA Unadjusted", "CSA Adjusted"))]
  
  hplot <- ggplot(dd_history, aes(x = event_time, y = est, color = model)) +
    geom_point(size = 2, position = position_dodge(width = 0.7)) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), 
                  width = 0, size = 1, position = position_dodge(width = 0.7)) +
    geom_hline(aes(yintercept = 0), linetype = 'dashed', size = 1) +
    scale_color_manual(values=c("#8856a7", "#2166ac"), 
                       breaks = c("CSA Unadjusted", "CSA Adjusted")) +
    facet_wrap(~outcome) +
    labs(title = "Effect of CFH on crime rates per 1k people",
         y = "Effect size", x = "Time relative to treatment", color = "Model") +
    theme_bw() +
    theme(legend.position = 'bottom',
          plot.title = element_text(hjust = 0.5),
          axis.title = element_text(family = "sans", size = 12),
          axis.title.y = element_text(angle = 0, vjust = 0.7),
          axis.text = element_text(family = "sans", size = 10),
          legend.text = element_text(family = 'sans', size = 10),
          strip.background = element_blank())
  
  ggsave(plot = hplot, filename = "./figs/model_results/csa/event_study.pdf",
         device = "pdf", height = 8.3, width = 11.7, units = "in")
  
}

