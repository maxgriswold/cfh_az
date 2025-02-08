# Prep covariates for AZ analysis
# Max Griswold
# 2/7/25

rm(list = ls())

library(data.table)
library(plyr)
library(dplyr)
library(tidycensus)

setwd("C:/users/griswold/Desktop/datasets/")

# Load concatennated files from Jacob Kaplan

#Kaplan, Jacob. Jacob Kaplan’s Concatenated Files: Uniform Crime Reporting Program Data: Offenses Known and Clearances by Arrest (Return A), 1960-2022. 
#Ann Arbor, MI: Inter-university Consortium for Political and Social Research [distributor], 2023-10-26. https://doi.org/10.3886/E100707V20

ucr <- 
ucr <- fread("ucr/ucr_1960_2020_offenses_known.csv")
