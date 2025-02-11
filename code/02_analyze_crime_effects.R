# Use augsynth to determine if crime-free housing led to reductions in crime
# Max Griswold
# 2/10/25

library(augsynth)
library(data.table)
library(ggplot2)
library(did)

df <- fread("./data/processed/df_crime_analysis.csv")

# Remove outlier locations, always-treated locations, and single pre-treatment locations

df <- df[outlier == F,]
df <- df[implementation_date > 2010|is.na(implementation_date),]

mod_1 <- multisynth(index_total ~ post_treat, location, year, data = df)
res_1 <- summary(mod_1)

mod_2 <- multisynth(burglary_total ~ post_treat, location, year, data = df)
res_2 <- summary(mod_2)

mod_3 <- multisynth(assault_total ~ post_treat, location, year, data = df)
res_3 <- summary(mod_3)

# Set separate synthetic controls for each cohort; see if results differ

mod_4 <- multisynth(index_total ~ post_treat, location, year, 
                    data = df, time_cohort = TRUE)
res_4 <- summary(mod_4)

mod_5 <- multisynth(burglary_total ~ post_treat, location, year, 
                    data = df, time_cohort = TRUE)
res_5 <- summary(mod_5)

mod_6 <- multisynth(assault_total ~ post_treat, location, year, 
                    data = df, time_cohort = TRUE)
res_6 <- summary(mod_6)

# Add covariates; using same covariates as those in the RAND report.

mod_7 <- multisynth(index_total ~ post_treat | pop_white + pop_black + pop_asian + pop_native_american + pop_latin_hispanic + median_income_adj, 
                    location, year, 
                    data = df, time_cohort = TRUE)
res_7 <- summary(mod_7)

mod_8 <- multisynth(burglary_total ~ post_treat | pop_white + pop_black + pop_asian + pop_native_american + pop_latin_hispanic + median_income_adj, location, year, 
                    data = df, time_cohort = TRUE)
res_8 <- summary(mod_8)

mod_9 <- multisynth(assault_total ~ post_treat | pop_white + pop_black + pop_asian + pop_native_american + pop_latin_hispanic + median_income_adj, location, year, 
                    data = df, time_cohort = TRUE)
res_9 <- summary(mod_9)

# Does using callaway and sant'anna change anything? 
df[is.na(implementation_date), implementation_date := 0]
df[, location := as.numeric(as.factor(location))]

csa_1 <- att_gt(yname = "index_total", 
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
