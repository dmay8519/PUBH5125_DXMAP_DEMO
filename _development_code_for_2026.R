#==============================================================================#
# Program:    
# Version:    1.0
# Author:     Dr Darren Mayne  
# Date:       
# Contact:    darren.mayne@health.nsw.gov.au  
# Purpose:    
# Notes:      
# Usage:      
#==============================================================================# 
# PARAMETERS: 
# -Name---------------- -Description-------------------------------------------# 
#                                             
#==============================================================================# 
# AMENDMENT HISTORY:
# -Ini- -Date-- -Id---- -Description-------------------------------------------# 
#                       
#==============================================================================# 
# This is public domain software. No guarantee as to suitability or accuracy is 
# given or implied. Users use this code entirely at their own risk.
#==============================================================================#

{ # Check for and install pacman package if missing
  
  if(base::require(pacman, quietly = TRUE) == FALSE){install.packages("pacman")}
  
  # list required packages
  
  reqPackages <- c("tidyverse",      # An opinionated collection of R packages designed for data science
                   "CARBayes",       # Spatial Generalised Linear Mixed Models for Areal Unit Data
                   "coda",           # A package for Output Analysis and Diagnostics for MCMC
                   "epitools",       # A package for epidemiological analysis
                   "ggmcmc",         # MCMC diagnostics with ggplot
                   "ggpubr",         # A package for producing publication-ready plots 
                   "ggrepel",        # A package for position non-overlapping text labels with 'ggplot2'
                   "GGally",         # A package to do pairplots
                   "knitr",          # A general-purpose tool for dynamic report generation in R
                   "nimble",         # A package for performing MCMC in R
                   "patchwork",      # A package to combine plots
                   "rmarkdown",      # A package for creating dynamic documents for R
                   "RColorBrewer",   # A package of colour palettes
                   "scales",         # Scale functions for visualization
                   "sf",             # Simple features for R
                   "spdep")          # A package to calculate neighbors)
  
  # Load (and install if missing) other required packages
  
  pacman::p_load(char = reqPackages, install = TRUE)
  
  # Check required packages are loaded
  
  pacman::p_loaded(char = reqPackages)
  
}

#------------------------------------------------------------------------------#
# Extract, transform and load (ETL) operations                              ----
#------------------------------------------------------------------------------#

# Set Census data and digital boundary URLs

abs_sa2 <- c("2021_GCP_SA2_for_NSW_short-header" = "https://www.abs.gov.au/census/find-census-data/datapacks/download/2021_GCP_SA2_for_NSW_short-header.zip",
             "SA2_2021_AUST_SHP_GDA2020" = "https://www.abs.gov.au/statistics/standards/australian-statistical-geography-standard-asgs-edition-3/jul2021-jun2026/access-and-downloads/digital-boundary-files/SA2_2021_AUST_SHP_GDA2020.zip")

# Download data

for (i in 1:base::length(abs_sa2)){
  
  utils::download.file(abs_sa2[i],
                       destfile = base::file.path(base::tempdir(),
                                                  base::paste0(base::names(abs_sa2[i]), ".zip"),
                                                  fsep = "\\"),
                       mode = "wb")
  
  utils::unzip(zipfile = base::file.path(base::tempdir(),
                               base::paste0(base::names(abs_sa2[i]), ".zip"),
                               fsep = "\\"),
               exdir = base::tempdir(),
               junkpaths = TRUE)
  
}

sa2_g19 <- utils::read.csv(file = base::paste(base::tempdir(),
                                              "2021Census_G19A_NSW_SA2.csv",
                                              sep = "\\")) %>%
  
  dplyr::inner_join(utils::read.csv(file = base::paste(base::tempdir(),
                                                      "2021Census_G19B_NSW_SA2.csv",
                                                      sep = "\\")),
                   by = "SA2_CODE_2021") %>%
  
  dplyr::inner_join(utils::read.csv(file = base::paste(base::tempdir(),
                                                      "2021Census_G19C_NSW_SA2.csv",
                                                      sep = "\\")),
                   by = "SA2_CODE_2021") %>%
  
  dplyr::rename_with(base::tolower) %>%
  
  tidyr::pivot_longer(cols = -c(sa2_code_2021),
                      names_to = "column",
                      values_to = "count") %>%
  
  dplyr::mutate(sex = base::factor(base::regmatches(column, base::regexpr("^[a-zA-z]", column)),
                                   levels = c("m", "f", "p"),
                                   labels = c("Males", "Females", "Persons")),
                age = base::regmatches(column, base::regexpr("([0-9]|(tot)){1, }.*$", column)),
                other = base::regmatches(column, base::regexpr("^[^a-zA-z]", column)))


utils::browseURL(base::tempdir())

column <- sa2_g19$column

base::regmatches(column, base::regexpr("[^a-zA-z].*", column))
