# Prep covariates for AZ analysis
# Max Griswold
# 2/7/25

rm(list = ls())

library(data.table)
library(plyr)
library(dplyr)
library(tidycensus)

setwd("C:/Users/griswold/Documents/GitHub/cfh_az/")

census_key <- "bdb4891f65609f274f701e92911b94365992028a"

###############
# Process ACS #
###############

extract_census <- function(y, geo = "place", us_state = "AZ"){
  
  # Grab tigris shapefiles w/ renter pop and total pop
  var_list <- data.table(variable_name = 
                           c("median_age", "median_age_male", "median_age_female",
                             "total_population", "pop_white", "pop_black", "pop_asian",
                             "pop_pacific_islander", "pop_native_american", "pop_other", "pop_immigrant",
                             "pop_male", "pop_female", "pop_latin_hispanic", "number_rental_units",
                             "number_housing_units", "number_vacant_units", "median_income", 
                             "agg_mortgage", "median_rent", "ratio_income_poverty_total",
                             "number_income_poverty_below_0.49", "number_income_poverty_0.5_to_0.99", 
                             "number_income_poverty_1_to_1.24", "number_income_poverty_1.25_to_1.49",
                             "pop_tenants", "pop_owners", "renter_household_total",
                             "renter_black", "renter_native_american", "renter_asian",
                             "renter_white", "renter_hispanic_latin", "renter_pacific_islander",
                             "renter_multiple_race", "renter_other", "renter_white_alone"),
                         variable_id   = 
                           c("B01002_001", "B01002_002", "B01002_003",
                             "B01003_001", "B02008_001", "B02009_001", "B02011_001",
                             "B02012_001", "B02010_001", "B02013_001", "B05012_003",
                             "B01001_002", "B01001_026", "B03001_003", "B25002_002",
                             "B25002_001", "B25002_003", "B19013_001",
                             "B25082_002", "B25064_001", "C17002_001",
                             "C17002_002", "C17002_003", 
                             "C17002_004", "C17002_005", 
                             "B25033_008", "B25003_002", "B25003_003",
                             "B25003B_003", "B25003C_003", "B25003D_003",
                             "B25003A_003", "B25003I_003", "B25003E_003",
                             "B25003G_003", "B25003F_003", "B25003H_003"))
  
  education_list <- data.table(variable_name = paste0("education=", rep(c(0, 2.5, 6.5, 7.5, 10, 11, 12, 12, 13, 14, 15, 15, 17, 19, 19, 21), 2)),
                               variable_id = sprintf("B15002_0%s", sprintf("%02d", c(seq(3, 18, 1), seq(20, 35, 1)))))
  
  var_list <- rbind(var_list, education_list)
  
  # Query API
  census_vars <- get_acs(geography = geo, year = y, variables = var_list$variable_id, 
                         state = us_state, geometry = F, key = census_key, survey = "acs5") %>%
                          setDT %>%
                          .[, .(GEOID, NAME, variable, estimate, moe)] %>%
                          setnames(., names(.), c("geoid", "location", "variable_id", "mean", "margin_of_error")) %>%
                          .[, `:=`(mean = as.numeric(mean),
                                   margin_of_error = as.numeric(margin_of_error))] %>%
                          join(., var_list, by = c("variable_id"), type = "left")

  #Process the education variables separately.
  education_vars <- census_vars[variable_name %like% "education", ]
  pop_vars       <- census_vars[variable_name %like% "pop",]
  renter_vars    <- census_vars[variable_name %like% "renter",]
  census_vars    <- census_vars[!(variable_name %like% "education" | variable_name %like% "pop" ), ]
  
  #Combine education variables into average years of schooling
  education_vars[, years := as.numeric(gsub(".*\\=", "", variable_name))]
  education_vars[, pop_years := sum(.SD$mean*.SD$years, na.rm = T), by = "location"]
  education_vars[, pop := sum(.SD$mean, na.rm = T), by = "location"]
  
  education_vars[, avg_edu := pop_years/pop]
  education_vars <- unique(education_vars[, .(location, avg_edu)])
  
  #Convert counts of renter housholds into percentages
  renter_vars <- dcast(renter_vars, location ~ variable_name, value.var = "mean")
  for (household in names(renter_vars)[names(renter_vars) %like% "renter" & !(names(renter_vars) %like% "total")]){
    
    renter_vars[, (household) := get(household)/renter_household_total]
    
  }
  
  #Convert counts of population into percentages
  pop_vars <- dcast(pop_vars, location ~ variable_name, value.var = "mean")
  for (populate in names(pop_vars)[names(pop_vars) %like% "pop" & !(names(pop_vars) %like% "total")]){
    
    pop_vars[, (populate) := get(populate)/total_population]
    
  }
  
  census_vars <- dcast(census_vars, geoid + location ~ variable_name, value.var = "mean")

  #Combine poverty percentage categories
  census_vars[, percent_poverty_150 := (number_income_poverty_below_0.49 + 
                                          number_income_poverty_0.5_to_0.99 + 
                                          number_income_poverty_1_to_1.24 + 
                                          number_income_poverty_1.25_to_1.49)/ratio_income_poverty_total]
  
  census_summary <- census_vars[, .(geoid, location, median_income, median_rent, 
                                    percent_poverty_150, number_rental_units)]
  
  census_summary <- join(census_summary, education_vars, by = "location", type = "left")
  census_summary <- join(census_summary, pop_vars, by = "location", type = "left")
  census_summary <- join(census_summary, renter_vars, by = "location", type = "left")
  
  census_summary[, year := y]

  return(census_summary)
  
}

years <- 2009:2023
census_vars <- setDT(ldply(years, extract_census))

# Adjust rent and income using CPI (2020 base year)

# U.S. Bureau of Labor Statistics, Consumer Price Index for All Urban Consumers: All Items in U.S. City Average [CPIAUCSL], 
# retrieved from FRED, Federal Reserve Bank of St. Louis; https://fred.stlouisfed.org/series/CPIAUCSL, February 7, 2025.

# U.S. Bureau of Labor Statistics, Consumer Price Index for All Urban Consumers: Rent of Primary Residence in U.S. City Average [CUUR0000SEHA], 
# retrieved from FRED, Federal Reserve Bank of St. Louis; https://fred.stlouisfed.org/series/CUUR0000SEHA, February 7, 2025.

# Using Rent of Primary residence due to inclusion of owner's equivalent rent within Shelter CPI

df_cpi_rent <- fread("./data/census/bls_cpi_urban_rent_average_non_seasonal.csv")
df_cpi_all  <- fread("./data/census/bls_cpi_urban_all_items_average_non_seasonal.csv")

# Recalculate CPI so that base year is 2023
base_2023_rent <- df_cpi_rent[year == 2023, cpi_rent]
df_cpi_rent[, cpi_rent := cpi_rent/base_2023_rent]

base_2023_all <- df_cpi_all[year == 2023, cpi_all]
df_cpi_all[, cpi_all := cpi_all/base_2023_all]

census_vars <- join(census_vars, df_cpi_rent, by = "year", type = 'left')
census_vars <- join(census_vars, df_cpi_all, by = "year", type = 'left')

census_vars[, median_rent_adj := median_rent*cpi_rent]
census_vars[, median_income_adj := median_income*cpi_all]

write.csv(census_vars, "./data/census/acs_az_place_2009_2023.csv", row.names = F)


###############
# Process UCR #
###############

# Load concatennated files from Jacob Kaplan

#Kaplan, Jacob. Jacob Kaplan’s Concatenated Files: Uniform Crime Reporting Program Data: Offenses Known and Clearances by Arrest (Return A), 1960-2022. 
#Ann Arbor, MI: Inter-university Consortium for Political and Social Research [distributor], 2023-10-26. https://doi.org/10.3886/E100707V20

ucr <- readRDS("./data/uniform_crime_reporting/offenses_known_yearly_1960_2023.rds")

setDT(ucr)
ucr <- ucr[state == "arizona" & year >= 2009]

# UCR distinguishes between two broad categories of crime: those committed
# against persons & those committed against property. Let's use these
# as the primary outcome variables, with their suggested aggregation in the
# IPCSR UCR documentation:

keep_vars <- c("ori", "year", "state", "crosswalk_agency_name",
               "number_of_months_missing", "geoid", 
               names(ucr)[names(ucr) %like% "actual_"])

# Reconstruct geoid from fips codes. Only hold onto the first six characters
# for GEOID, which correspond to CDP (essentially, we're upscaling more
# granular geographies to be at the CDP-level). Make sure to include
# leading zeros so that fips place code is always six characters

ucr[, fips_place_code := formatC(fips_place_code, width = 5,
                                 format = "d", flag = "0")]

ucr[, geoid := paste0(fips_state_code, fips_place_code)]
ucr[, geoid := sub("^(\\d{6})", "\\1", geoid)]

ucr <- ucr[, keep_vars, with = F]

# Remove locations missing a geoid or agency_name
ucr <- ucr[(!((geoid %like% "NA")|crosswalk_agency_name == ""))]

# Melt data by crime_type, then calculate totals for each crime type by 
# census place. 
ucr <- melt(ucr, id.vars = c("ori", "year", "state", "number_of_months_missing",
                             "geoid", "crosswalk_agency_name"), 
            variable.name = "crime_type", value.name = "crime_count")

# Only hold onto specific crime types
ucr[, crime_type := gsub("actual_", "", crime_type)]

# Collapse counts by year, CDP, crime-type. Also aggregate number of LEAs informing
# counts, along with cumulative months missing:
ucr[, number_of_reporting_agencies := .N, by = c("crime_type", "geoid", "year")]
ucr[, number_of_months_missing := sum(.SD$number_of_months_missing, na.rm = T), by = c("crime_type", "geoid", "year")]

ucr[, crime_count := sum(.SD$crime_count, na.rm = T), by = c("crime_type", "geoid", "year")]
ucr <- setDT(unique(ucr[, .(geoid, year, crime_type, crime_count, 
                            number_of_months_missing, number_of_reporting_agencies)]))

write.csv(ucr, "offenses_known_az_2009_2023.csv", row.names = F)

######################
# Merge census & ucr #
######################

# Get restricted set 
pop_vars <- census_vars[, .(geoid, location, total_population, year)]

# Reshape UCR long so each crime category is a separate column
ucr <- dcast(ucr, geoid + year + number_of_months_missing + number_of_reporting_agencies ~ crime_type, value.var = "crime_count")

test <- setDT(join(pop_vars, ucr, by = c("geoid", "year"), type = "right"))
test <- test[!is.na(all_crimes)]
