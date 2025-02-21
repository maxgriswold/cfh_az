# Prep covariates for AZ analysis
# Max Griswold
# 2/7/25

rm(list = ls())

library(data.table)
library(plyr)
library(dplyr)
library(tidycensus)
library(tidygeocoder)
library(sf)

sf::sf_use_s2(FALSE)
options(tigris_use_cache = TRUE)

setwd("C:/Users/griswold/Documents/GitHub/cfh_az/")

census_key <- "bdb4891f65609f274f701e92911b94365992028a"

overwrite <- F
plots     <- F
geocode   <- F

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

if (overwrite){

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
  
  # Clean up location names to better merge with treatment data
  census_vars[, location := gsub(", Arizona", "", location)]
  census_vars[, location := gsub(" CDP| town| city", "", location)]

  write.csv(census_vars, "./data/census/acs_az_place_2009_2023.csv", row.names = F)
}else{
  census_vars <- fread("./data/census/acs_az_place_2009_2023.csv")
}

###############
# Process UCR #
###############

# Load concatennated files from Jacob Kaplan

#Kaplan, Jacob. Jacob Kaplan’s Concatenated Files: Uniform Crime Reporting Program Data: Offenses Known and Clearances by Arrest (Return A), 1960-2022. 
#Ann Arbor, MI: Inter-university Consortium for Political and Social Research [distributor], 2023-10-26. https://doi.org/10.3886/E100707V20

if (overwrite){

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
  
  write.csv(ucr, "./data/uniform_crime_reporting/offenses_known_az_2009_2023.csv", row.names = F)
  
}else{
  ucr <- fread("./data/uniform_crime_reporting/offenses_known_az_2009_2023.csv")
}

######################
# Merge census & ucr #
######################

# Get restricted set 
pop_vars <- census_vars[, .(geoid, location, total_population, year)]

ucr <- setDT(join(ucr, pop_vars, by = c("geoid", "year"), type = "left"))

# Drop all-crimes variable and locations which are not one of the 91 AZ cities
ucr <- ucr[crime_type != "all_crimes"]
ucr <- ucr[!is.na(location)]

# Convert crime counts to rate per 1k
ucr[, crime_rate_1k := (crime_count/total_population)*1000]

# Break out into a few key categories. Aggregating myself, rather than using
# FBI provided aggregate categories, since the index total miss some key crim
# types. For example, the assault_index excludes assaults committed with a gun

ucr[, crime_subtype := crime_type]
ucr[, crime_type := NA]

# Remove total categories:
ucr[crime_subtype %like% c("burg"), crime_type := 'burglary_total']
ucr[crime_subtype %like% c("robbery"), crime_type := 'robbery_total']
ucr[crime_subtype %like% c("assault"), crime_type := 'assault_total']
ucr[crime_subtype %like% c("mtr"), crime_type := 'motor_theft_total']
ucr[crime_subtype %like% c('ind'), crime_type := "index_total"]

ucr <- ucr[!is.na(crime_type)]

# Investigate data to determine if any cities should be outliered in crime analysis
if (plots){
  crime_type_by_year <- function(type, city, dd){
    
    dd <- dd[location == city, ]
    
    if (type != "pop" & type != "missingness"){
      dd <-  dd[crime_type == type,]
    }
    
    if (type != "pop" & type != "missingness"){
      
      plot <- ggplot(dd, aes(x = year, y = crime_rate_1k, color = crime_subtype)) +
        geom_line(linewidth = 1, alpha = 0.9) +
        labs(title = paste0(type, " rates per 10k in ", city),
             x = "Year",
             y = "Rate per 10k") +
        theme_bw() +
        theme(strip.background = element_blank())
      
    }else if (type == "pop"){
      
      dd <- unique(dd[, .(year, total_population)])
      
      plot <- ggplot(dd, aes(x = year, y = total_population)) +
        geom_line(linewidth = 1) +
        labs(title = paste0("Population estimates in ", city),
             x = "Year",
             y = "Population") +
        theme_bw() +
        theme(strip.background = element_blank())
      
    }else{
      
      dd <- unique(dd[, .(year, number_of_months_missing, number_of_reporting_agencies)])
      
      plot <- ggplot(dd, aes(x = year, y = number_of_months_missing)) +
        geom_col() +
        labs(title = paste0("Number of months not reported across ", 
                            unique(dd$number_of_reporting_agencies),
                            " agencies"),
             x = "Year",
             y = "Number of months") +
        theme_bw() +
        theme(strip.background = element_blank())
      
    }
    
    return(plot)
    
  }
  
  args <- expand.grid("type" = c(unique(ucr$crime_type), c("pop", "missingness")), 
                      "city" = unique(ucr$location))
  
  crimes_year <- mlply(args, crime_type_by_year, dd = ucr)
  crimes_year <- marrangeGrob(crimes_year, nrow = 3, ncol = 2)
  
  ggsave("./figs/crimes_by_city_year.pdf", crimes_year, height = 21, width = 29.7, units = "cm")
  
}

# Investigating the plots, I'm noticing the following irregularities:

# St. Johns did not report data before 2012
# Benson did not report several years before 2016
# Bisbee did not report 2014-2016
# Douglas did not report 2013 - 2015
# Huachuca did not report 2013 - 2016
# Sierra Vista missing 2013
# Tombstone missing 2012 - 2016
# Sedona - One agency did not report data before 2016; time trends look reasonable otherwise. 2021 appears to be missing
# Globe is missing 2012
# Hayden is missing data before 2016
# Winkelman is super small - outlier this!
# Safford is missing data between 2010 - 2014
# Thatcher has one agency missing data before 2016 - trends appear relatively reasonable (besides 2021)
# Clifton has implausible trends - does not appear to be reporting data fter 2016
# Outlier Duncan
# Outlier Gila Bend - they do not report any crime!
# Mesa is missing data from one agency before 2016 but otherwise looks reasonable
# Peoria is missing data from two agencies before 2016 but otherwise looks reasonable
# Tolleson is missing data in 2016
# Outlier Youngtown
# Outlier Guadalupe
# Outlier Colorado City
# Outlier Taylor
# South Tucson is missing data in 2016 (and 2017)
# Sells data looks implausible
# Superior is missing data before 2016
# Outlier Patagonia
# San Luis is missing data before 2011
# Wellton is missing data in 2012 - 2013
# Parker is missing data from 1/3 agencies in 2013 - 2016
# Quartzsite is missing data in 2012

# Implausible large values:
# Window Rock, San Carlos, White River, Peach Springs, and Sacaton are reporting extremely large crime rates
# for their population size.

# Based on the above, I will be outliering these locations from the crime analysis. This comprises
# a total population of 36k/5.6 million in AZ CDPs.

outlier_locs <- c("Hayden", "Clifton", "Duncan", "Youngtown", "Guadalupe", "Colorado City", "Taylor", "Superior",
                  "Sells", "Patagonia", "Tombstone", "Benson", "Sacaton", "San Carlos", "Window Rock", "Whiteriver",
                  "Peach Springs", "Bisbee", "Gila Bend", "Safford", "South Tucson", "Winkelman")

# Remove total categories and "simple" categories (reason for aggregating myself rather
# than using index - assault totals exclude assaults with a gun)

ucr <- ucr[!(crime_subtype %like% "assault_total"|crime_subtype %like% "burg_total"|crime_subtype %like% "robbery_total"|
                       crime_subtype %like% "mtr_veh_theft_total")]

ucr[, crime_rate_1k := sum(.SD$crime_rate_1k), by = c("geoid", "year", "crime_type")]

ucr <- setDT(unique(ucr[, .(geoid, location, year, number_of_months_missing, number_of_reporting_agencies, crime_type, crime_rate_1k)]))

ucr <- dcast(ucr, geoid + location + year + number_of_months_missing + number_of_reporting_agencies ~ crime_type,
                 value.var = "crime_rate_1k")

ucr <- ucr[, outlier := ifelse(location %in% outlier_locs, T, F)]

# Add on census information and treatment data for final analysis dataset:
df_analysis <- setDT(join(ucr, census_vars, by = c("geoid", "location", "year")))

# We're missing one city which isn't reporting any information to UCR/NIBRS:
# Queen Creek
# so add on census statistics and outlier location for crime analysis

df_analysis <- rbind(df_analysis, census_vars[location == "Queen Creek"], fill = T)
df_analysis[location == "Queen Creek", outlier := T]

treated_locs <- fread("./data/crime_free_housing/cfh_treated_sites.csv")

df_analysis <- setDT(join(df_analysis, treated_locs, by = c("location"), type = "left"))
df_analysis[, treated := ifelse(is.na(treated), F, T)]

# Recode treatment status so augsynth runs correctly. 
df_analysis[, post_treat := ifelse(year >= implementation_date & treated == T, 1, 0)]

# For now, remove Bullhead as crime-free since I cannot find evidence it exists:
df_analysis[location == "Bullhead City", `:=`(treated = F, post_treat = 0)]

write.csv(df_analysis, "./data/processed/df_crime_analysis.csv", row.names = F)

if (plots){
  
  df_prepost <- df_analysis[, .(year, assault_total, burglary_total, index_total, treated)]
  df_prepost <- melt(df_prepost, id.vars = c("year", "treated"), 
                     variable.name = "crime_type", value.name = "crime_rate_1k")
  
  df_prepost <- df_prepost[!is.na(crime_rate_1k)]
  df_prepost[, crime_rate_1k := mean(.SD$crime_rate_1k), by = c("crime_type", "treated", "year")]
  
  df_prepost <- unique(df_prepost)
  
  ggplot(df_prepost, aes(x = year, y = crime_rate_1k, color = treated)) +
    geom_line() +
    facet_wrap(~crime_type) +
    theme_bw()
  
}

############################
# Process Eviction Filings #
############################

df_evict <- fread("./data/eviction_filings/maricopa_county/evictions_maricopa_county_2024.csv")

df_evict <- df_evict[, c("File Date", "Case Number", "SubCategory", "Case Status",
                         "Claim Amount", "Defendant1 Address", "Plaintiff Name"), with = F]

setnames(df_evict, names(df_evict), c("date", "case_id", "category", "status", "rent_owed", "address", "plaintiff"))

# Hold onto only adjudicated cases in 2024. Large number of cases are sealed (16k)
# which indicates either the defendant won the case or both parties agreed to seal
# the case.
df_evict <- df_evict[status != "Sealed",]

# Hold onto possession cases (this is the overwhelming number of cases in Maricopa county:
# 65k/70k are possession. ~5k are rent-related)

df_evict <- df_evict[category == "Poss Prop",]

# Are cases uniquely identified by an ID?
# dim(df_evict)[1] == length(unique(df_evict$case_id))

batch_geocodes <- function(d, batches = 2e3, add = ""){
  
  #Create splitting column
  d  <- d[, batch := ceiling(.I/batches)]
  ds <- split(d, by = "batch", keep.by = F)
  
  processed <- llply(ds, code_coordinates, added_text = add)
  processed <- rbindlist(processed)
  
  return(processed)
  
}

# Iterative function which applies geocoders to addresses. For each address,
# try a geocoder. If it doesn't work, try the next one.

code_coordinates <- function(df, added_text = ""){
  
  #Geocoders are super temperamental. Given this, we'll need to process geocodes
  #in batches. We'll split datasets into sets of 2k. For each batch, we'll try
  #out a cascade of open-source geocoders, starting w/ Census, then arc-gis, 
  #then finally open-street-maps. If all else fails, we can try Google at a later point (I have
  #an API key but it's set up to bill a different PTN)
  
  #Lots of unit tests in this function; I was running into  issues
  #from certain geocoders not producing new results. There's likely a better
  #way to do the below; too much rewritten code!
  
  dd <- df %>%
    geocode(address, method = "census") %>%
    setDT() %>%
    .[, method := "census"]
  
  #Add on additional text to address to create additional redundancies, which
  #tends to improve reliability of geocodes
  dd[, address := paste(address, added_text)]
  
  res <- dd[!is.na(lat)]
  dd  <- dd[is.na(lat),]
  
  #Remove lat/long column for uncoded addresses since keeping this column
  #in the dataset is creating issues for arcgis method:
  dd <- dd[, `:=`(lat = NULL, long = NULL)]
  
  if (dim(dd)[1] >= 1){
    dd <- dd %>%
      geocode(address, method = "arcgis") %>%
      setDT() %>%
      .[, method := "arcgis"]
    
    if (dim(dd[!is.na(lat)])[1] >= 1){
      res <-setDT(rbind(res, dd[!is.na(lat)]))
    }
    dd  <- dd[is.na(lat), ]
    dd <- dd[, `:=`(lat = NULL, long = NULL)]
  }
  
  if (dim(dd)[1] >= 1){
    dd <- dd %>%
      geocode(address, method = "osm") %>%
      setDT() %>%
      .[, method := "osm"]
    
    if (dim(dd[!is.na(lat)])[1] >= 1){
      res <- setDT(rbind(res, dd[!is.na(lat)]))
    }
  }
  
  print(paste0("Could not code ", dim(dd)[1], " addresses."))
  print(res)
  return(res)
  
}

# Also load CFH property locations within cities of Maricopa county.
# Geocode these locations.

if (geocode){
  
  keep_cities <- c("avondale", "chandler", "mesa")
  
  property_loc_files <- "./data/crime_free_housing/property_locations/"
  property_locs <- paste0(property_loc_files, list.files(property_loc_files))
  property_locs <- property_locs[grepl(paste0(keep_cities, collapse = "|"), property_locs)]
  
  code_properties <- function(city){
    
    city_name <- gsub(paste0(property_loc_files, "|.csv"), "", city)
    
    dd <- fread(city)
    
    # Add on property name, if it exists, to help out geocoder:
    if ("property_name" %in% names(dd)){
      dd[, address := paste0(property_name, ", ", address) ]
    }
    
    dd <- batch_geocodes(dd, add = paste0(", ", city_name, ", AZ"))
    dd[, location := city_name]
    
    return(dd)
    
  }
  
  df_properties <- ldply(property_locs, code_properties)
  write.csv(df_properties, "./data/crime_free_housing/cfh_property_locations.csv")
  
}else{
  
  df_properties <- fread("./data/crime_free_housing/cfh_property_locations.csv")
  
}


if (geocode){
  
  df_evict <- batch_geocodes(df_evict)
  
  write.csv(df_evict, "./data/eviction_filings/evictions_maricopa_county_2024_geocoded.csv", row.names = F)
  
  # Convert eviction points into counts at block-group level:
  df_evict <- geo_in_block(df_evict, tiger_block)
  df_evict[, evict_count := .N, by = "block_geoid"]
  df_evict <- unique(df_evict[, .(block_geoid, evict_count)])
  
  write.csv(df_evict, "./data/eviction_filings/evictions_maricopa_county_2024_block_groups.csv", row.names = F)
  
}else{
  
  df_evict <- fread("./data/eviction_filings/evictions_maricopa_county_2024_block_groups.csv")
  df_evict[, block_geoid := as.numeric(block_geoid)]
}

# Get TIGER Files for both census places and blocks. 
# We'll use TIGER places to crop blocks into cities, then
# merge on eviction filings and property locations, if they 
# intersect a block.

# Variable we're querying doesn't matter - we're doing this to get geometry
# only.

tiger_state <- get_acs(geography = "state", year = 2023, variables = "B01003_001", 
                       state = "AZ", geometry = T, key = census_key, survey = "acs5") %>%
                setDT %>%
                .[, .(GEOID, NAME, geometry)] %>%
                setnames(., c("NAME", "GEOID"), c("block_name", "block_geoid")) %>%
                .[, block_geoid := as.numeric(block_geoid)] %>%
                st_as_sf()

tiger_county <- get_acs(geography = "county", year = 2023, variables = "B01003_001", 
                       state = "AZ", geometry = T, key = census_key, survey = "acs5") %>%
                        setDT %>%
                        .[, .(GEOID, NAME, geometry)] %>%
                        setnames(., c("NAME", "GEOID"), c("block_name", "block_geoid")) %>%
                        .[, block_geoid := as.numeric(block_geoid)] %>%
                        st_as_sf()

tiger_place <- get_acs(geography = "place", year = 2023, variables = "B01003_001", 
                       state = "AZ", geometry = T, key = census_key, survey = "acs5") %>%
                        setDT %>%
                        .[, .(GEOID, NAME, geometry)] %>%
                        setnames(., c("NAME", "GEOID"), c("place_name", "place_geoid")) %>%
                        .[, place_geoid := as.numeric(place_geoid)] %>%
                        st_as_sf()

tiger_block <- get_acs(geography = "block group", year = 2023, variables = "B01003_001", 
                       state = "AZ", geometry = T, key = census_key, survey = "acs5") %>%
                        setDT %>%
                        .[, .(GEOID, NAME, geometry)] %>%
                        setnames(., c("NAME", "GEOID"), c("block_name", "block_geoid")) %>%
                        .[, block_geoid := as.numeric(block_geoid)] %>%
                        st_as_sf()

st_write(tiger_state, "./data/shapefiles/tiger_state.geojson", append = F, 
         delete_dsn = T, delete_layer = T)

st_write(tiger_county, "./data/shapefiles/tiger_county.geojson", append = F, 
         delete_dsn = T, delete_layer = T)

tiger_place$location <- tolower(gsub(" CDP| town| city|, Arizona", "", tiger_place$place_name))
st_write(tiger_place, "./data/shapefiles/tiger_place.geojson", append = F, 
         delete_dsn = T, delete_layer = T)

# Hold onto TIGER places within study sites (e.g., within Maricopa county)

tiger_place <- tiger_place[tiger_place$location %in% keep_cities,]
setorder(tiger_place, place_name)

# Using TIGER shapefile, determine all blocks that touch a place polygon
# Then, for each block, calculate the percentage of area within a census 
# place. 

block_in_place <- function(place, blocks){
  
  overlap               <- st_intersection(place, blocks)
  overlap$area_overlap  <- as.numeric(st_area(overlap))
  
  # Intersections can occur at the boundary between two the city & neighboring
  # city blocks. Drop these
  overlap <- overlap[overlap$area_overlap > 0,]
  
  setDT(overlap)
  overlap <- overlap[, .(block_geoid, area_overlap)]
  
  # Get original area of blocks and merge with intersected blocks
  old_blocks <- tiger_block[tiger_block$block_geoid %in% overlap$block_geoid,]
  old_blocks$area_original  <- as.numeric(st_area(old_blocks))
  
  setDT(old_blocks)
  old_blocks <- old_blocks[, .(block_geoid, area_original, geometry)]
  
  overlap <- join(overlap, old_blocks, by = "block_geoid", type = "left")
  overlap[, overlap_percent := area_overlap/area_original, ]
  overlap[, overlap_80 := ifelse(overlap_percent >= 0.8, 1, 0)]
  overlap[, overlap_50 := ifelse(overlap_percent >= 0.5, 1, 0)]
  
  return(overlap)
  
}

# For a given spatial point, determine which block it belongs to.
geo_in_block <- function(d, geo){
  
  d <- st_as_sf(d, coords = c("long", "lat"))
  d <- st_set_crs(d, st_crs(tiger_block))
  
  d <- st_join(d, geo, join = st_within, left = T)
  setDT(d)
  
  d[, block_geoid := as.numeric(block_geoid)]
  d[, .(address, block_geoid)]
  
  return(d)
  
}

# Below code looks a bit goofy but the idea here is to transform the single
# sf collection of cities into a list where each item is now a sf corresponding
# to a single row in the original collection. I am doing this so I can then
# loop over the shapefile of each city to determine which census blocks are within
# them.

city_shape  <- lapply(seq_len(nrow(tiger_place)), function(i) tiger_place[i, , drop = F])
city_blocks <- lapply(city_shape, block_in_place, blocks = tiger_block)

df_census_blocks <- extract_census(2023, geo = "block group")
df_census_blocks <- df_census_blocks[location %like% "Maricopa",]
df_census_blocks[, block_geoid := as.numeric(geoid)]

# Get all files in alignment and combine assorted information across files
names(city_shape) <- keep_cities
names(city_blocks) <- keep_cities

sfd_analysis_prep <- function(city_name){
  
  dd <- city_blocks[[city_name]]
  
  # Add on census indicators, eviction information, and cfho property
  # locations:
  dd <- join(dd, df_census_blocks, by = "block_geoid", type = 'left')
  dd <- join(dd, df_evict, by = "block_geoid", type = "left")
  
  dd_props <- df_properties[location == city_name,]
  dd_props <- geo_in_block(dd_props, tiger_block)
  dd_props[, cfh_num := .N, by = "block_geoid"]
  dd_props <- unique(dd_props[, .(block_geoid, cfh_num)])
  
  dd <- join(dd, dd_props, by = "block_geoid", type = "left")
  
  # If blocks are missing eviction or cfh property counts, then
  # we know they should be zero:
  dd[is.na(evict_count), evict_count := 0]
  dd[is.na(cfh_num), cfh_num := 0]
  
  dd[, evict_rate := 100*(evict_count/number_rental_units)]
  dd[, median_income_10k := median_income/1e4]
  
  dd[, cfh_any := ifelse(cfh_num > 0, 1, 0)]  
  
  # Do not include blocks without rental units as a comparators
  dd <- dd[number_rental_units > 0,]
  
  dd[, location := city_name]
  
  dd <- dd[, .(location, block_geoid, evict_count, evict_rate,
               cfh_num, cfh_any, number_rental_units,
               renter_household_total, total_population, median_income_10k,
               renter_white_alone, renter_black, renter_asian, renter_native_american, 
               renter_hispanic_latin, percent_poverty_150, geometry)]
  
  dd <- st_as_sf(dd) %>%
          st_transform(4326)
  
  sfd_save <- paste0("./data/processed/sfd/", city_name, "_prepped.geojson")
  st_write(dd, sfd_save, append = F, delete_dsn = T, delete_layer = T)
  
  return(dd)
  
}

shapefiles <- lapply(keep_cities, sfd_analysis_prep)

# Sanity check - how do the shapes overlap?
# mesa_city   <- city_shape[[3]]
# mesa_blocks <- st_as_sf(city_blocks[[3]]) %>% st_set_crs(., st_crs(mesa_city))
# mesa_props  <- st_as_sf(df_properties[df_properties$location == 'mesa',], coords = c("long", "lat")) %>% st_set_crs(., st_crs(mesa_city))
# ggplot() + geom_sf(data = mesa_city, color = 'red', fill = "red", alpha = 0.3) +
#            geom_sf(data = mesa_blocks, fill = "transparent") +
#            geom_sf(data = mesa_props)

# Looks correct to me.
