#==============================================================================#
# Program:    disease_mapping_demonstration_code.R
# Version:    1.0
# Author:     Dr Darren Mayne  
# Date:       
# Contact:    darren.mayne@sydney.edu.au
# Purpose:    This R script demonstrates how to fit the Besag, York, and Mollié
#             (1991) model for disease mapping and ecological regression
#             applications. It uses a classic data set describing male lip
#             cancer incidence in Scottish administration districts from
#             1975-1980. The data are available from the GeoDa Centre GuitHub
#             repository (see link below), which also provides original and
#             supplementary references for the data.
# Notes:      Session information
#
#             R version 4.1.3 (2022-03-10)
#             Platform: x86_64-w64-mingw32/x64 (64-bit)
#             Running under: Windows 10 x64 (build 18363)
# 
#             attached base packages:
#               stats
#               graphics
#               grDevices
#               utils
#               datasets
#               methods
#               base
# 
#             other attached packages:
#              spdep_1.3-3
#              spData_2.2.2
#              sf_1.0-16
#              scales_1.3.0
#              RColorBrewer_1.1-3
#              patchwork_1.2.0
#              nimble_1.1.0
#              ggrepel_0.9.5
#              ggpubr_0.6.0
#              ggmcmc_1.5.1.1
#              epitools_0.5-10.1
#              coda_0.19-4.1
#              lubridate_1.9.2
#              forcats_1.0.0
#              stringr_1.5.0
#              dplyr_1.1.0       
#              purrr_1.0.1
#              readr_2.1.4
#              tidyr_1.3.0
#              tibble_3.2.0
#              ggplot2_3.5.1
#              tidyverse_2.0.0 
#
#             Check package versions above if executing code results in errors.
#             ggplot2_3.5.0 or above must be installed to produce maps
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

base::setwd(base::tempdir())

### Session set up -------------------------------------------------------------


# User pacman package to load packages needed for demonstration

{ # Check for and install pacman package if missing
  
  if(base::require(pacman, quietly = TRUE) == FALSE){install.packages("pacman")}
  
  # list required packages
  
  reqPackages <- c("tidyverse",      # A package for data manipulation and plotting
                   "coda",           # A package for Output Analysis and Diagnostics for MCMC
                   "epitools",       # A package for epidemiological analysis
                   "ggmcmc",         # MCMC diagnostics with ggplot
                   "ggpubr",         # A package for producing publication-ready plots 
                   "ggrepel",        # A package for position non-overlapping text labels with 'ggplot2'
                   "GGally",         # A package to do pairplots
                   "nimble",         # A package for performing MCMC in R
                   "patchwork",      # A package to combine plots
                   "RColorBrewer",   # A package of colour palettes
                   "scales",         # Scale functions for visualization
                   "sf",             # Simple features for R
                   "spdep")          # A package to calculate neighbors)
  
  # Load (and install if missing) other required packages
  
  pacman::p_load(char = reqPackages, install = TRUE)
  
  # Check required packages are loaded
  
  pacman::p_loaded(char = reqPackages)
}

# Get session information
#
# utils::sessionInfo()

## Set session constants -------------------------------------------------------

# Create a vector of district names to exclude when labeling maps

exc_labels <- c("Eastwood",
                "Kilmarnock",
                "Monklands",
                "Motherwell",
                "Renfrew",
                "East Kilbride",
                "Cumbernauld",
                "Kirkcaldy",
                "West Lothian",
                "Inverclyde",
                "Strathkelvin",
                "Bearsden",
                "East Lothian",
                "Midlothian")

# Visualization colour palette

col_pal <- c("low" = RColorBrewer::brewer.pal(11, "RdYlBu")[11],
             "mid" = RColorBrewer::brewer.pal(11, "RdYlBu")[6],
             "high" = RColorBrewer::brewer.pal(11, "RdYlBu")[1])

#------------------------------------------------------------------------------#
# Data preparation                                                             #
#------------------------------------------------------------------------------#

### Extract data ---------------------------------------------------------------

## Download Scottish lip cancer data for demonstration  ------------------------

# Available from the the GeoDa Center as a zipped archive from:
# https://geodacenter.github.io/data-and-lab/scotlip/. However, we can download
# it directly within R using the following code chunk.

data.tmp <- base::file.path(tempdir(), "scotlip.zip")
data.url <- "https://geodacenter.github.io/data-and-lab/data/scotlip.zip"
utils::download.file(data.url, data.tmp, mode = "wb")

# View the contents of the zip file

utils::unzip(data.tmp, list = TRUE) %>%
  base::as.data.frame() %>%
  arrange(Name)

## Extract the scotlip.gpkg geopackage data file (row 79) ----------------------

# See https://en.wikipedia.org/wiki/GeoPackage for a description of the 
# geopackage data format. 

# Extract the scotlip.gpkg geopackage the extracted ZIP file

utils::unzip(data.tmp,
             files = "scotlip/scotlip.gpkg",
             exdir = "Data",
             junkpaths = TRUE) # Only use the base name of the stored file path when extracting

# List GeoPack layers (there is only 1 called "scotlip")

base::data.frame(layers = sf::st_layers("Data/scotlip.gpkg")$name)

# Read scotlip.gpkg layer, which includes data and geometry.

# Note: I have set the coordinate reference system (CRS) to NA because the CRS
#       supplied in the geopack (WGS 80) is incorrect.

src_sf <- sf::st_read(dsn = "Data/scotlip.gpkg",
                      layer = "scotlip",
                      crs = NA) # Remove the wrong coordinate reference system (CRS) assigned in the geopackage

# Attach the correct CRS for the layer: OSGB36 / British National Grid

sf::st_crs(src_sf) <- 27700 # see https://epsg.io/27700 for CRS details

## Use gpk_sf layer to create a Scotland boundary layer ----------------------

bnd_sf <- sf::st_union(src_sf, by_feature = FALSE) %>%
  sf::st_as_sf() %>%
  dplyr::mutate(area = "Scotland") %>%
  dplyr::select(area, geom = x)

## View data embedded in geopack -----------------------------------------------

utils::head(src_sf, n = 10)

### Transform data set ---------------------------------------------------------

# Create an analytic sf for data analysis

lip_sf <- src_sf %>%
  
  # Rename columns to lowercase
  
  dplyr::rename_with(base::tolower) %>%
  
  # Limit to selected columns
  
  dplyr::select(dist_code = district, dist_name = name, pop, obs = cancer, exp = cexp, aff) %>%
  
  # Correct typographical errors in district names
  
  dplyr::mutate(dplyr::across(dist_name, ~ base::gsub("EastKilbride", "East Kilbride",
                                                      base::gsub("EastLothian", "East Lothian",
                                                                 base::gsub("NEFife", "NE Fife",
                                                                            base::gsub("WesternIsles", "Western Isles",
                                                                                       base::gsub("WestLothian", "West Lothian",.)))))),
                one_over_exp = 1 / exp) %>%
  
  # Calculate standardised incidence ratios (SIR) with exact Poisson 95% 
  # confidence intervals using using pois.exact function from epitools.
  
  dplyr::bind_cols(base::data.frame(epitools::pois.exact(pull(., obs), pull(., exp)) %>%
                                      dplyr::select(obs_sir_est = rate,
                                                    obs_sir_l95 = lower,
                                                    obs_sir_u95  = upper))) %>%
  
  # Calculate standard and relative errors for SIR
  
  dplyr::mutate(obs_sir_se = base::sqrt(obs) / exp,
                obs_sir_rse = obs_sir_se / obs_sir_est) %>% # or 1 / base::sqrt(obs)

  # Add variable labels
  
  labelled::set_variable_labels(dist_code = "District code",
                                dist_name = "District name",
                                pop = "Male population years at risk (1975-1980)",
                                obs = "Observed male lip cancer cases (1975-1980)",
                                exp = "Expected male lip cancer cases (1975-1980)",
                                aff = "Percent of male population working in Agriculture, Fishing or Forestry",
                                one_over_exp = "Proportional relative standard error based on EXP",
                                obs_sir_est = "Observed SIR estimate",
                                obs_sir_l95 = "Observed SIR lower 95% CI",
                                obs_sir_u95 = "Observed SIR upper 95% CI",
                                obs_sir_se = "Observed SIR standard error",
                                obs_sir_rse = "Observed SIR relative standard error",
                                geom = "Geometry") %>%
  
  # Move geometry column to end of data set
  
  dplyr::select(dplyr::everything(), -geom, geom)

## View analytic data set ------------------------------------------------------

View(lip_sf)

#------------------------------------------------------------------------------#
# Initial map of male lip cancer risk in Scotland, 1975-1980                   #
#------------------------------------------------------------------------------#

### Map observed SIR -----------------------------------------------------------

# Find the absolute maximum SIR on the log scale

max_log_sir <- base::as.integer(base::max(base::abs(base::log(lip_sf[lip_sf$obs != 0, ]$obs_sir_est))) / 1) + 1

# Produce the SIR map

map_obs_sir <- ggplot2::ggplot() +
  
  # Add districts with SIR > 0
  
  ggplot2::geom_sf(data = lip_sf %>% dplyr::filter(obs_sir_est > 0),
                   mapping = ggplot2::aes(fill = obs_sir_est),
                   col = "grey50") +
  
  # Add districts with SIR = 0
  
  ggplot2::geom_sf(data = lip_sf %>% dplyr::filter(obs_sir_est == 0),
                   mapping = ggplot2::aes(alpha = base::factor("SIR = 0", ordered = TRUE)),
                   fill = "grey90",
                   colour = "grey50",
                   size = 0.1) +
  
  # Add Scotland boundary
  
  ggplot2::geom_sf(data = bnd_sf,
                   fill = NA,
                   colour = "Black",
                   size = 0.25) +
  
  # Add fill aesthetic gradient mapping
  
  ggplot2::scale_fill_gradient2("Observed SIR",
                                low = col_pal["low"], # blue
                                mid = col_pal["mid"], # yellow
                                high = col_pal["high"], # darkblue
                                midpoint = 0,
                                na.value = NA,
                                limits = c(base::exp(-max_log_sir), base::exp(max_log_sir)),
                                breaks = base::exp(base::seq(-max_log_sir, max_log_sir, 0.5)),
                                labels = base::format(base::round(exp(seq(-max_log_sir, max_log_sir, 0.5)), 2), nsmall = 2),
                                trans = "log") +
  
  ggplot2::scale_alpha_manual(values = c("SIR = 0" = 1)) +
  
  # Add code labels to identify districts
  
  ggrepel::geom_text_repel(data = lip_sf %>% dplyr::filter(! dist_name %in% exc_labels),
                           inherit.aes = FALSE,
                           mapping = ggplot2::aes(label = dist_name,
                                                  geometry = geom),
                           stat = "sf_coordinates",
                           force = 0,
                           color = "white",
                           bg.color = "grey30",
                           bg.r = 0.1,
                           size = 9/.pt) +
  
  # Pretty up the legend
  
  ggplot2::guides(fill = ggplot2::guide_colourbar(frame.colour = "black", 
                                                  frame.linewidth = 0.25,
                                                  ticks.colour = "black",
                                                  ticks.linewidth = 0.25,
                                                  order = 1),
                  alpha = ggplot2::guide_legend(title = NULL, order = 2)) +
  
  # Pretty up the map
  
  ggplot2::theme_bw() + 
  
  ggplot2::theme(axis.title = ggplot2::element_blank(),
                 axis.text = ggplot2::element_blank(),
                 axis.ticks = ggplot2::element_blank(),
                 
                 legend.position = "inside",
                 legend.margin = ggplot2::margin(-3, 0, 0, 0),
                 legend.position.inside = c(0.15, 0.88),
                 
                 panel.grid = ggplot2::element_blank(),
                 
                 text = ggplot2::element_text(size = 12))

### Map relative standard errors for the observed SIR --------------------------

map_obs_rse <- ggplot2::ggplot(data = lip_sf) +
  
  # Add relative standard errors for districts with SIR > 0
  
  # Add districts with SIR > 0
  
  ggplot2::geom_sf(data = lip_sf %>% dplyr::filter(obs_sir_est > 0),
                   mapping = ggplot2::aes(fill = obs_sir_rse),
                   col = "grey50") +
  
  # Add districts with SIR = 0
  
  ggplot2::geom_sf(data = lip_sf %>% dplyr::filter(obs_sir_est == 0),
                   mapping = ggplot2::aes(alpha = base::factor("SIR = 0", ordered = TRUE)),
                   fill = "grey90",
                   colour = "grey50",
                   size = 0.1) +
  
  # Add code labels to identify districts
  
  ggrepel::geom_text_repel(data = lip_sf %>% dplyr::filter(! dist_name %in% exc_labels),
                           inherit.aes = FALSE,
                           mapping = ggplot2::aes(label = dist_name,
                                                  geometry = geom),
                           stat = "sf_coordinates",
                           force = 0,
                           color = "white",
                           bg.color = "grey30",
                           bg.r = 0.1,
                           size = 9/.pt) +
  
  # Add fill aesthetic gradient mapping
  
  ggplot2::scale_fill_gradient2("Relative standard error",
                                low = col_pal["low"], # blue
                                mid = col_pal["mid"], # yellow
                                high = col_pal["high"], # darkblue
                                midpoint = 0.5,
                                na.value = NA,
                                limits = c(0, 1),
                                breaks = base::seq(0, 1, 0.1),
                                labels = base::seq(0, 1, 0.1)) +
  
  ggplot2::scale_alpha_manual(values = c("SIR = 0" = 1)) +
  
  ggplot2::guides(fill = ggplot2::guide_colourbar(frame.colour = "black", 
                                                  frame.linewidth = 0.25,
                                                  ticks.colour = "black",
                                                  ticks.linewidth = 0.25,
                                                  order = 1),
                  alpha = ggplot2::guide_legend(title = NULL, order = 2)) +
  
  ggplot2::theme_bw() + 
  
  ggplot2::theme(axis.title = ggplot2::element_blank(),
                 axis.text = ggplot2::element_blank(),
                 axis.ticks = ggplot2::element_blank(),
                 
                 legend.position = "inside",
                 legend.margin = ggplot2::margin(-3, 0, 0, 0),
                 # legend.position.inside = c(0.15, 0.88),
                 legend.position.inside = c(0.20, 0.88),
                 legend.text = ggplot2::element_text(size = 10),
                 
                 panel.grid = ggplot2::element_blank(),
                 
                 text = ggplot2::element_text(size = 12))

## Output paired observed SIR and expected standard errors ---------------------

map_sir_rse <- (map_obs_sir | map_obs_rse) +
  patchwork::plot_annotation(title = "Observed standardised incidence ratios (SIR) and relative standard errors (RSE) for male lip cancer in Scotland (1975-1980)", 
                             theme = theme(plot.title = element_text(hjust = 0.5)))

base::print(map_sir_rse)

# Save map - not run
# 
# ggplot2::ggsave(map_sir_rse,
#                 filename = "Figures/Observed standardised incidence ratios (SIR) and relative standard errors (RSE) for male lip cancer in Scotland (1975-1980).png",
#                 device = "png",
#                 height = 99.3 * 4,
#                 width = 228.6 * 4,
#                 units = "mm",
#                 bg = "white",
#                 dpi = 300)

#------------------------------------------------------------------------------#
# Create Queen adjacency (contiguity) list for each district                   #
#------------------------------------------------------------------------------#

### Create Queen adjacency (contiguity) list for each district -----------------

# The spdep::poly2nb function constructs a list of neighbours from a polygon object

adj_nb_list <- spdep::poly2nb(lip_sf, 
                              queen = TRUE, 
                              row.names = row.names(lip_sf))

# View adjacency list structure

utils::head(adj_nb_list)

### Output adjacency matrix in WinBUGS format for use in BYM models ------------

# The spdep::nb2WB function can be used to make WinBUGS neighbours list 

adj_nb_wgt <- spdep::nb2WB(nb = adj_nb_list)

# View WinBUIGS adjacency list

base::print(adj_nb_wgt)

# The *adj* vector defines the district indexes for the neighbours of the ith
# analysis unit; the *weights* vector indicates the weight assigned to each
# neighbour when calculating the conditional autoregression spatial term
# in the BYM model; and *num* indicates how many of the neighbours in *adj*
# border the ith analysis unit. For example, the first element of *num* indicates
# that the district of Skye-Lochalsh (row 1 in *lip_sf*) has three neighbours, 
# which *adj* indicates correspond to row indexes 5, 9 and 19 (Ross-Cromarty,
# Lochaber and Inverness districts) in *lip_sf*, and are assigned weights of 1, 
# 1 and 1 when calculating the CAR spatial term for Skye-Lochalsh.

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#                                                                              #
# Smooth SIR estimates using BYM disease mapping model                         #
#                                                                              #
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

### Assemble objects needed to run a BYM disease model -------------------------

## Obtain the number of districts ----------------------------------------------

N <- dim(lip_sf)[1]

## Format SIR data for NIMBLE in a list ----------------------------------------

bym_data = base::list(obs = lip_sf$obs)               # Observed lip cancer cases

## List constants required by NIMBLE ti fit BYM model --------------------------

bym_consts <-base::list(N = N,                          # Number of districts 
                        
                        exp = lip_sf$exp,               # Expected lip cancer cases
                        
                        # Adjacency matrix values used to specify CAR distribution
                        
                        L = length(adj_nb_wgt$weights),  # Total number of neighbours        
                        adj = adj_nb_wgt$adj,            # Vector of neighbours for each ditrict                            
                        num = adj_nb_wgt$num,            # Vector of neighbour counts for each district
                        weights = adj_nb_wgt$weights)    # Weight for each neighbour (all 1 for a queen matrix

## Write the model specification -----------------------------------------------

# NOTE: Don't use R namespaces as models are run by the nimble C++ compiler

bym_code <- nimble::nimbleCode(
  {
    for (i in 1:N){
      
      # Poisson likelihood for observed counts
      
      obs[i] ~ dpois(lambda[i]) 
      log(lambda[i]) <- alpha + s[i] + u[i] + log(exp[i])
      
      # Prior for the area-specific unstructured random effect
      
      u[i] ~ dnorm(0, tau = tau.u)    
      
      # Smoothed area-specific SIR
      
      smoothed.sir[i] <- exp(alpha + s[i] + u[i])
      
      # Overall residual area-specific cancer risk above or below the Scottish average
      
      residual.sir[i] <- exp(s[i] + u[i])
      
      # Residual area-specific cancer risk due to unobserved and spatially structured factors (s)
      
      residual.sir.s[i] <- exp(s[i])
      
      # Residual area-specific lip cancer risk due to unobserved and spatially unstructured factors (u)
      
      residual.sir.u[i] <- exp(u[i])
      
      # Posterior probabilities for residual risk - is the SIR for the ith district higher or lower than 1
      
      residual.sir.gt1[i] <- 1 - step(1 - residual.sir[i]) # Posterior probability that SIR > 1
      residual.sir.lt1[i] <- 1 - step(residual.sir[i] - 1) # Posterior probability that SIR < 1                         
    } 
    
    # Prior for the area-specific spatially structured random effect (ICAR)
    
    s[1:N] ~ dcar_normal(adj[1:L],
                         weights[1:L],
                         num[1:N],
                         tau.s,
                         zero_mean = 1)
    
    # Vague uniform prior for intercept
    
    alpha ~ dflat()                                       
    
    # Overall SIR for the study region (lip cancer risk common to all areas)
    
    overall.sir <- exp(alpha)                                
    
    # Hyperprior distribution on inverse variance (precision) for the unstructured variance component (u)
    
    tau.u ~ dgamma(1, 0.01)                            
    
    # Variance of u (unstructured variance component)
    
    sigma2.u <- 1/tau.u                                   
    
    # Hyperprior distribution on inverse variance (precision) for the spatially structured variance component (s)
    
    tau.s ~ dgamma(1, 0.01)                              
    
    # variance of s (spatial variance component)
    
    sigma2.s <- 1/tau.s
    
  }
)

## Specify initial values for nodes as a data list -----------------------------

# alpha = intercept
# u     = prior for unstructured random effect
# tau.u = hyperprior for the precision of u
# s     = prior for structured random effect
# tau.s = hyperprior for the precision of s

bym_inits <- list(list(alpha = 0.01,                              # MCMC chain 1 
                       u = base::rep(0.01, times = N), 
                       tau.u = 10,
                       s = base::rep(0.01, times = N),
                       tau.s = 10),
                  list(alpha = 0.5,                               # MCMC chain 2
                       u = base::rep(-0.01, times = N),
                       tau.u = 1,
                       s = base::rep(-0.01, times = N),
                       tau.s = 1))

## Select parameters (nodes) to monitor ----------------------------------------

bym_params <- c("alpha",
                "u",
                "s",
                "sigma2.s",
                "sigma2.u",
                "overall.sir",
                "residual.sir",
                "residual.sir.s",
                "residual.sir.u",
                "smoothed.sir",
                "residual.sir.gt1",
                "residual.sir.lt1")

## Specify the MCMC sampler settings -------------------------------------------

n_chains <- 2             # Number of MCMC chains to run  
n_burnin <- 50000         # Number of iterations to use as model burn-in
n_iter <- n_burnin * 2    # Number of iterations per chain
n_thin <- 10              # Thinning interval for chain, i.e., keep every nth sample

# Burn-in samples per chain (n = 2) - used to converge model

(n_burnin / n_thin)                # 5000 samples used to burn the model in

# Inference samples per chain (n = 2)

((n_iter - n_burnin) / n_thin)     # 5000 samples from the posterior distribution summarised for inference

### Fit the model --------------------------------------------------------------

# The nimble showCompilerOutput and nimbleMCMC options are set to TRUE for
# demonstration purposes only. You can exclude these options if you prefer

nimbleOptions(showCompilerOutput = TRUE)
bym_samples <- nimble::nimbleMCMC(code = bym_code,
                                  data = bym_data,
                                  constants = bym_consts, 
                                  inits = bym_inits,
                                  monitors = bym_params,
                                  niter = n_iter,
                                  nburnin = n_burnin,
                                  thin = n_thin, 
                                  nchains = n_chains, 
                                  setSeed = 9,              # Seed for randome number generator - important for reproducibility  
                                  progressBar = TRUE,
                                  samplesAsCodaMCMC = TRUE, 
                                  summary = TRUE, 
                                  WAIC = TRUE)
nimbleOptions(showCompilerOutput = FALSE)

# See https://stats.stackexchange.com/q/304958 for discussion about WAIC warning

## View raw list output returned my the MCMC sampler ---------------------------

# Sampled values from the posterior distribution (n = n_iter - n_burin = 5,000)

View(bym_samples$samples$chain1)     
View(bym_samples$samples$chain2)

# Summarised sample values from the posterior distribution (mean, median, std and quantiles 2.5 and 97.5)

View(bym_samples$summary$chain1)
View(bym_samples$summary$chain2)
View(bym_samples$summary$all.chains)

### Check model convergence using the Gelman-Rubin convergence diagnostic ------

# The Gelman-Rubin compares the within chain variation to the between chain
# variation for each parameter. If the chains have converged on the stationary
# distribution, then it should be impossible to distinguish between the chains.
# It is a ratio measures, and values < 1.1 are considered evidence of chain
# convergence.

gr.diag <- coda::gelman.diag(bym_samples$samples, multivariate = FALSE)

# List all parameters that are not posterior probabilities

gr.diag.inc <- stringr::str_subset(base::dimnames(gr.diag$psrf)[[1]], "gt1|lt1", negate = TRUE)

## Are all values < 1.1? If yes, this is evidence of convergence ---------------

base::all(gr.diag$psrf[gr.diag.inc,"Point est."] < 1.1) # Don't check posterior probabilities - these are threshold nodes and do not need to converge

## If no, then Which parameters have a value >= 1.1? ---------------------------

base::which(gr.diag$psrf[gr.diag.inc,"Point est."] > 1.1)

### Check trace, density and autocorrelation plots for convergence -------------

# Import MCMC samples into a ggs that can be used with ggs_* graphical functions

bym_ggmcmc <- ggmcmc::ggs(bym_samples$samples)## Check the overall SIR -------------------------------------------------------

# Generate trace, density and autocorrelation plots

trace_overall.sir <- bym_ggmcmc %>% 
  dplyr::filter(Parameter == "overall.sir") %>% 
  ggmcmc::ggs_traceplot() +
  ggplot2::theme_bw() +
  ggplot2::theme(text = ggplot2::element_text(size = 12))

density_overall.sir <- bym_ggmcmc %>% 
  dplyr::filter(Parameter == "overall.sir") %>% 
  ggmcmc::ggs_density() +
  ggplot2::theme_bw() +
  ggplot2::theme(text = ggplot2::element_text(size = 12))

acf_overall.sir <- bym_ggmcmc %>%
  dplyr::filter(Parameter == "overall.sir") %>% 
  ggmcmc::ggs_autocorrelation() +
  ggplot2::theme_bw() +
  ggplot2::theme(text = element_text(size = 12))

# Plot visual convergence diagnostics

base::print((trace_overall.sir | density_overall.sir) / acf_overall.sir) +
  patchwork::plot_annotation(title = "Trace, density and autocorrelation plots for overall SIR",
                             theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)))

## Check the smoothed SIR for Perth-Kinross (29) -------------------------------

# This district has the most neighbours

# Generate trace, density and autocorrelation plots

trace_smoothed.sir_29 <- bym_ggmcmc %>% 
  dplyr::filter(Parameter == "smoothed.sir[29]") %>% 
  ggmcmc::ggs_traceplot() +
  ggplot2::theme_bw() +
  ggplot2::theme(text = ggplot2::element_text(size = 12))

density_smoothed.sir_29 <- bym_ggmcmc %>% 
  dplyr::filter(Parameter == "smoothed.sir[29]") %>% 
  ggmcmc::ggs_density() +
  ggplot2::theme_bw() +
  ggplot2::theme(text = ggplot2::element_text(size = 12))

acf_smoothed.sir_29 <- bym_ggmcmc %>%
  filter(Parameter == "smoothed.sir[29]") %>% 
  ggmcmc::ggs_autocorrelation() +
  ggplot2::theme_bw() +
  ggplot2::theme(text = ggplot2::element_text(size = 12))

# Plot visual convergence diagnostics

base::print((trace_smoothed.sir_29 | density_smoothed.sir_29) / acf_smoothed.sir_29) +
  patchwork::plot_annotation(title = "Trace, density and autocorrelation plots for Perth-Kinross smoothed SIR",
                             theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)))

## Check the smoothed SIR for Shetland (8) -------------------------------------

# This district has no neighbours

# Generate trace, density and autocorrelation plots

trace_smoothed.sir_08 <- bym_ggmcmc %>% 
  dplyr::filter(Parameter == "smoothed.sir[8]") %>% 
  ggmcmc::ggs_traceplot() +
  ggplot2::theme_bw() +
  ggplot2::theme(text = ggplot2::element_text(size = 12))

density_smoothed.sir_08 <- bym_ggmcmc %>% 
  dplyr::filter(Parameter == "smoothed.sir[8]") %>% 
  ggmcmc::ggs_density() +
  ggplot2::theme_bw() +
  ggplot2::theme(text = ggplot2::element_text(size = 12))

acf_smoothed.sir_08 <- bym_ggmcmc %>%
  filter(Parameter == "smoothed.sir[8]") %>% 
  ggmcmc::ggs_autocorrelation() +
  ggplot2::theme_bw() +
  ggplot2::theme(text = ggplot2::element_text(size = 12))

# Plot visual convergence diagnostics

base::print((trace_smoothed.sir_08 | density_smoothed.sir_08) / acf_smoothed.sir_08) +
  patchwork::plot_annotation(title = "Trace, density and autocorrelation plots for Shetland smoothed SIR",
                             theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)))

## Check the spatially variance component --------------------------------------

# Generate trace, density and autocorrelation plots

trace_sigma2.s <- bym_ggmcmc %>% 
  dplyr::filter(Parameter == "sigma2.s") %>% 
  ggmcmc::ggs_traceplot() +
  ggplot2::theme_bw() +
  ggplot2::theme(text = ggplot2::element_text(size = 12))

density_sigma2.s <- bym_ggmcmc %>% 
  dplyr::filter(Parameter == "sigma2.s") %>% 
  ggmcmc::ggs_density() +
  ggplot2::theme_bw() +
  ggplot2::theme(text = ggplot2::element_text(size = 12))

acf_sigma2.s <- bym_ggmcmc %>%
  filter(Parameter == "sigma2.s") %>% 
  ggmcmc::ggs_autocorrelation() +
  ggplot2::theme_bw() +
  ggplot2::theme(text = ggplot2::element_text(size = 12))

# Plot visual convergence diagnostics

base::print((trace_sigma2.s | density_sigma2.s)/acf_sigma2.s) +
  patchwork::plot_annotation(title = "Trace, density and autocorrelation plots for spatially correlated variance component (s)",
                             theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)))

## Check the uncorrelated variance component -----------------------------------

# Generate trace, density and autocorrelation plots

trace_sigma2.u <- bym_ggmcmc %>% 
  dplyr::filter(Parameter == "sigma2.u") %>% 
  ggmcmc::ggs_traceplot() +
  ggplot2::theme_bw() +
  ggplot2::theme(text = ggplot2::element_text(size = 12))

density_sigma2.u <- bym_ggmcmc %>% 
  dplyr::filter(Parameter == "sigma2.u") %>% 
  ggmcmc::ggs_density() +
  ggplot2::theme_bw() +
  ggplot2::theme(text = ggplot2::element_text(size = 12))

acf_sigma2.u <- bym_ggmcmc %>%
  filter(Parameter == "sigma2.u") %>% 
  ggmcmc::ggs_autocorrelation() +
  ggplot2::theme_bw() +
  ggplot2::theme(text = ggplot2::element_text(size = 12))

base::print((trace_sigma2.u | density_sigma2.u) /acf_sigma2.u) +
  patchwork::plot_annotation(title = "Trace, density and autocorrelation plots for uncorrelated variance component (u)",
                             theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)))

## Summary of convergence checking ---------------------------------------------

# Gelman-Rubin diagnostic provides evidence for convergence for all parameters;
# however, trace plots indicate poor convergence for the unstructured variance
# component (sigma2.u). If this was not a demonstration analysis, I would tune
# the MCMC settings to improve the convergence of the unstructured variance
# parameter; however, we will just proceed to summarising the posterior samples

### Summarise across MCMC samples to obtained BYM model estimates --------------

## Extract alpha, sigma2.s and sigma2.u ----------------------------------------

bym_ests <- bym_samples$summary$all.chains[c("alpha",
                                             "sigma2.s",
                                             "sigma2.u"),
                                           c("Median",
                                             "95%CI_low",
                                             "95%CI_upp")] %>%
  
  # Convert matrix to data frame
  
  base::as.data.frame() %>%
  
  # Rename columns
  
  stats::setNames(c("median", "lower_95_ci", "upper_95_ci")) %>%
  
  # Make the data frame row names an explicit variable called *parameter*
  
  tibble::rownames_to_column(var = "parameter") %>%
  
  # Exponentiate the parameter and 95% CI estimates for alpha to obtain overall SIR effect measure
  
  dplyr::mutate(dplyr::across(c(median, lower_95_ci, upper_95_ci), ~ base::ifelse(parameter == "alpha", base::exp(.), NA), .names = "{.col}_exp")) %>%
  
  # Format estimates as text variables with three decimal places
  
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ base::ifelse(base::is.na(.), "", base::format(base::round(., 3), nsmall = 3))))

### Print model estimates and WIAC ---------------------------------------------

base::print(list("Model estimate" = bym_ests,
                 "WAIC" = bym_samples$WAIC$WAIC))

### Summarise across MCMC samples to obtain smoothed SIR estimates -------------

# Create data set of smoothed and residual risk of lip cancer

bym_sf <- lip_sf %>%
  
  # Smoothed SIR estimate
  
  dplyr::bind_cols(bym_samples$summary$all.chains[base::paste0("smoothed.sir[", 1:N, "]"), c("Median", "95%CI_low", "95%CI_upp")] %>%
                     base::as.data.frame() %>%
                     stats::setNames(c("bym_smoothed_sir_est", "bym_smoothed_sir_l95", "bym_smoothed_sir_uci"))) %>%
  
  # Residual SIR estimate
  
  dplyr::bind_cols(bym_samples$summary$all.chains[base::paste0("residual.sir[", 1:N, "]"), c("Median", "95%CI_low", "95%CI_upp")] %>%
                     base::as.data.frame() %>%
                     stats::setNames(c("bym_residual_sir_est", "bym_residual_sir_l95", "bym_residual_sir_uci"))) %>%
  
  # Residual s SIR estimate
  
  dplyr::bind_cols(bym_samples$summary$all.chains[base::paste0("residual.sir.s[", 1:N, "]"), c("Median", "95%CI_low", "95%CI_upp")] %>%
                     base::as.data.frame() %>%
                     stats::setNames(c("bym_residual_sir_s_est", "bym_residual_sir_s_l95", "bym_residual_sir_s_uci"))) %>%
  
  # Residual u SIR estimate
  
  dplyr::bind_cols(bym_samples$summary$all.chains[base::paste0("residual.sir.u[", 1:N, "]"), c("Median", "95%CI_low", "95%CI_upp")] %>%
                     base::as.data.frame() %>%
                     stats::setNames(c("bym_residual_sir_u_est", "bym_residual_sir_u_l95", "bym_residual_sir_u_uci"))) %>%
  
  # Posterior probabilities for residual risk > 1
  
  dplyr::bind_cols(bym_samples$summary$all.chains[base::paste0("residual.sir.gt1[", 1:N, "]"), "Mean"] %>%
                     base::as.data.frame() %>%
                     stats::setNames("bym_residual_sir_gt1")) %>%
  
  # Posterior probabilities for residual risk < 1
  
  dplyr::bind_cols(bym_samples$summary$all.chains[base::paste0("residual.sir.lt1[", 1:N, "]"), "Mean"] %>%
                     base::as.data.frame() %>%
                     stats::setNames("bym_residual_sir_lt1")) %>%
  
  # Annotate exceedence probabilities
  
  dplyr::mutate(bym_residual_sir_cat = dplyr::case_when(bym_residual_sir_gt1 > 0.995 ~ "++",
                                                        bym_residual_sir_gt1 > 0.975 ~ "+",
                                                        bym_residual_sir_lt1 > 0.995 ~ "--",
                                                        bym_residual_sir_lt1 > 0.975 ~ "-",
                                                        .default = ""))

### Map BYM smoothed SIR for districts -----------------------------------------

map_bym_residual_sir <- ggplot() +
  
  # Add SIR estimates
  
  geom_sf(data = bym_sf,
          aes(fill = bym_residual_sir_est),
          col = "grey50",
          linewidth = 0.1) +
  
  # Add SIR posterior probabilities
  
  ggplot2::geom_point(data = bym_sf %>% dplyr::filter(bym_residual_sir_cat != ""),
                      mapping = ggplot2::aes(shape = bym_residual_sir_cat,
                                             geometry = geom),
                      fill = "black",
                      stat = "sf_coordinates") +
  
  # Add Scotland boundary
  
  ggplot2::geom_sf(data = bnd_sf,
                   fill = NA,
                   colour = "Black",
                   linewidth = 0.25) +
  
  # Add fill aesthetic gradient mapping
  
  ggplot2::scale_fill_gradient2("Smoothed SIR",
                                low = col_pal["low"], # blue
                                mid = col_pal["mid"], # yellow
                                high = col_pal["high"], # darkblue
                                midpoint = 0,
                                na.value = NA,
                                limits = c(base::exp(-max_log_sir), base::exp(max_log_sir)),
                                breaks = base::exp(base::seq(-max_log_sir, max_log_sir, 0.5)),
                                labels = base::format(base::round(exp(seq(-max_log_sir, max_log_sir, 0.5)), 2), nsmall = 2),
                                trans = "log") +
  
  # Add shape aesthic mapping
  
  ggplot2::scale_shape_manual("Posterior probabilities",
                              values = c("++" = 24,
                                         "+" = 2,
                                         "-" = 6,
                                         "--" = 25),
                              breaks = c("++", "+", "-", "--"),
                              labels = c("++" = "SIR > 1, p > 97.5%",
                                         "+" = "SIR > 1, p > 95.0%",
                                         "-" = "SIR < 1, p > 95.0%",
                                         "--" = "SIR < 1, p > 97.5%")) +
  
  # Format legends
  
  guides(fill = guide_colourbar(frame.colour = "black", 
                                frame.linewidth = 0.25,
                                ticks.colour = "black",
                                ticks.linewidth = 0.25,
                                order = 1),
         shape = guide_legend(order = 2)) +
  
  # Format map
  
  theme_bw() + 
  
  theme(axis.title = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        
        legend.position = "inside",
        legend.margin = ggplot2::margin(0, 0, 0, 0),
        legend.position.inside = c(0.15, 0.88),
        
        panel.grid = ggplot2::element_blank(),
        
        text = ggplot2::element_text(size = 12))

## Output observed and BYM smoothed SIR for lip cancer risk --------------------

map_obs_bym <- (map_obs_sir | map_bym_residual_sir) +
  patchwork::plot_annotation(title = "Observed and BYM smoothed SIR for male lip cancer in Scotland (1975-1980)", 
                             theme = theme(plot.title = element_text(hjust = 0.5)))

base::print(map_obs_bym)

# Save map - not run
# 
# ggplot2::ggsave(map_obs_bym,
#                 filename = "Figures/Observed and BYM smoothed SIR for male lip cancer in Scotland (1975-1980).png",
#                 device = "png",
#                 height = 99.3 * 4,
#                 width = 228.6 * 4,
#                 units = "mm",
#                 bg = "white",
#                 dpi = 300)
                
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#                                                                              #
# Ecological regression application of BYM disease mapping model               #
#                                                                              #
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

# QUESTION: Is district-level lip cancer risk associated with percentage of male
#           population employed in agriculture, fishing and forestry (AFF)
#           industries? AFF serves as a proxy measure for sun exposure, which is
#           a risk factor for lip cancer.

### Map percentages of male populations in AFF industries ----------------------

map_aff <- ggplot2::ggplot(data = lip_sf) +
  
  # Add relative standard errors for districts with SIR > 0
  
  # Add districts with SIR > 0
  
  ggplot2::geom_sf(data = lip_sf,
                   mapping = ggplot2::aes(fill = aff),
                   col = "grey50") +
  
  # Add code labels to identify districts
  
  ggrepel::geom_text_repel(data = lip_sf %>% dplyr::filter(! dist_name %in% exc_labels),
                           inherit.aes = FALSE,
                           mapping = ggplot2::aes(label = dist_name,
                                                  geometry = geom),
                           stat = "sf_coordinates",
                           force = 0,
                           color = "white",
                           bg.color = "grey30",
                           bg.r = 0.1,
                           size = 9/.pt) +
  
  # Add fill aesthetic gradient mapping
  
  ggplot2::scale_fill_gradient2("Percent in AFF",
                                low = col_pal["low"], # blue
                                mid = col_pal["mid"], # yellow
                                high = col_pal["high"], # darkblue
                                midpoint = 12.5,
                                na.value = NA,
                                limits = c(0, 25),
                                breaks = base::seq(0, 25, 5),
                                labels = base::seq(0, 25, 5)) +
  
  ggplot2::guides(fill = ggplot2::guide_colourbar(frame.colour = "black", 
                                                  frame.linewidth = 0.25,
                                                  ticks.colour = "black",
                                                  ticks.linewidth = 0.25,
                                                  order = 1)) +
  
  ggplot2::theme_bw() + 
  
  ggplot2::theme(axis.title = ggplot2::element_blank(),
                 axis.text = ggplot2::element_blank(),
                 axis.ticks = ggplot2::element_blank(),
                 
                 legend.position = "inside",
                 legend.margin = ggplot2::margin(-3, 0, 0, 0),
                 # legend.position.inside = c(0.15, 0.88),
                 legend.position.inside = c(0.20, 0.88),
                 legend.text = ggplot2::element_text(size = 10),
                 
                 panel.grid = ggplot2::element_blank(),
                 
                 text = ggplot2::element_text(size = 12))

## Compare percentage of males in AFF versus smoothed lip cancer risk ----------

map_bym_aff <- (map_aff | map_bym_residual_sir) +
              patchwork::plot_annotation(title = "Percentage of male population in agriculture, fishing or forestry (AFF) versus smoothed risk for male lip cancer in Scotland (1975-1980)", 
                                         theme = theme(plot.title = element_text(hjust = 0.5)))

base::print(map_bym_aff)

# Save map - not run
# 
# ggplot2::ggsave(map_bym_aff,
#                 filename = "Figures/Percentage of male population in AFF versus smoothed risk for male lip cancer in Scotland (1975-1980).png",
#                 device = "png",
#                 height = 99.3 * 4,
#                 width = 228.6 * 4,
#                 units = "mm",
#                 bg = "white",
#                 dpi = 300)

### Assemble objects needed to run a BYM ecological regression -----------------

## Obtain the number of districts ----------------------------------------------

N <- dim(lip_sf)[1]

# Format SIR and AFF data for NIMBLE in a list ---------------------------------

eco_data = list(obs = lip_sf$obs,                         # Observed lip cancer cases
                aff = lip_sf$aff)                         # Percent male population in agriculture, fishing or forestry

# display data

base::print(eco_data)

## Format constants data for NIMBLE in a list ----------------------------------

eco_consts <-list(N = N,                                  # Number of districts 
                  
                  exp = lip_sf$exp,                       # Expected lip cancer cases
                  
                  # Adjacency matrix values used to specify CAR distribution
                  
                  L = length(adj_nb_wgt$weights),         # Total number of neighbours        
                  adj = adj_nb_wgt$adj,                   # Vector of neighbours for each ditrict                            
                  num = adj_nb_wgt$num,                   # Vector of neighbour counts for each district
                  weights = adj_nb_wgt$weights)           # Weight for each neighbour (all 1 for a queen matrix)

## Write the model specification -----------------------------------------------

# NOTE: Don't use R namespaces as models are run by the nimble C++ compiler - R is just the interface

eco_code <- nimble::nimbleCode(
  {
    for (i in 1:N){
      
      # Poisson likelihood for observed counts
      
      obs[i] ~ dpois(lambda[i]) 
      log(lambda[i]) <- log(exp[i]) + alpha + beta * aff[i] + s[i] + u[i]
      
      # Area-specific unstructured random effect
      
      u[i] ~ dnorm(0, tau = tau.u)    
      
      # Area-specific SIR
      
      smoothed.sir[i] <- exp(alpha + beta * aff[i] + s[i] + u[i])
      
      # Area-specific residual SIR
      
      residual.sir[i] <- exp(s[i] + u[i])
      
      # Area-specific residual SIR due to unobserved and spatially structured factors
      
      residual.sir.s[i] <- exp(s[i])
      
      # Area-specific residual SIR due to unobserved and unstructured factors
      
      residual.sir.u[i] <- exp(u[i])
      
      # Posterior probabilities
      
      residual.sir.gt1[i] <- 1 - step(1 - residual.sir[i]) # Posterior probability that SIR > 1
      residual.sir.lt1[i] <- 1 - step(residual.sir[i] - 1) # Posterior probability that SIR < 1
    } 
    
    # Prior for the area-specific spatially structured random effect (ICAR)
    
    s[1:N] ~ dcar_normal(adj[1:L],
                         weights[1:L],
                         num[1:N],
                         tau.s,
                         zero_mean = 1)
    
    # Vague uniform prior for intercept
    
    alpha ~ dflat() 
    
    # Overall SIR for the study region (lip cancer risk common to all areas)
    
    overall.sir <- exp(alpha)   
    
    # Normal distribution prior for beta (AFF effect measure)
    
    beta ~ dnorm(0, 1)
    
    # Risk ratio for a 1% increase in the percentage of males in AFF occupations
    
    aff.sir <- exp(beta)
    
    # Hyperprior distribution on inverse variance (precision) for the unstructured variance component (u)
    
    tau.u ~ dgamma(1, 0.01)                            
    
    # Variance of u (unstructured variance component)
    
    sigma2.u <- 1/tau.u                                   
    
    # Hyperprior distribution on inverse variance (precision) for the spatially structured variance component (s)
    
    tau.s ~ dgamma(1, 0.01)                              
    
    # variance of s (spatial variance component)
    
    sigma2.s <- 1/tau.s 
    
  }
  
)                                       

## Specify initial values for nodes --------------------------------------------

# alpha = intercept
# beta  = effect estimate for AFF
# u     = prior for unstructured random effect
# tau.u = hyper prior for the precision of u
# s     = prior for structured random effect
# tau.s = hyper prior for the precision of s

eco_inits <- list(list(alpha = 0.01,                             # MCMC chain 1 
                       beta = 1,
                       tau.u = 10,
                       u = base::rep(0.01, times = N), 
                       tau.s = 10,
                       s = base::rep(0.01, times = N)),
                  list(alpha = 0.5,                              # MCMC chain 2
                       beta = -1,
                       tau.u = 1,
                       u = base::rep(-0.01, times = N), 
                       tau.s = 1,
                       s = base::rep(-0.01, times = N)))

## Select parameters (nodes) to monitor ----------------------------------------

eco_params <- c("alpha",
                "beta",
                "u",
                "s",
                "sigma2.s",
                "sigma2.u",
                "overall.sir",
                "aff.sir",
                "residual.sir",
                "residual.sir.s",
                "residual.sir.u",
                "smoothed.sir",
                "residual.sir.gt1",
                "residual.sir.lt1")

## Specify the MCMC sampler settings -------------------------------------------

n_chains <- 2             # Number of MCMC chains to run 
n_burnin <- 50000         # Number of iterations to use as model burn-in
n_iter <- n_burnin * 2    # Number of iterations samples 
n_thin <- 10              # Thinning interval for chain, i.e., keep every nth sample

### Fit the model --------------------------------------------------------------

nimbleOptions(showCompilerOutput = TRUE)
eco_samples <- nimble::nimbleMCMC(code = eco_code,
                                  data = eco_data,
                                  constants = eco_consts, 
                                  inits = eco_inits,
                                  monitors = eco_params,
                                  niter = n_iter,
                                  nburnin = n_burnin,
                                  thin = n_thin, 
                                  nchains = n_chains, 
                                  setSeed = 9,              # Seed for randome number generator - important for reproducibility  
                                  progressBar = TRUE,
                                  samplesAsCodaMCMC = TRUE, 
                                  summary = TRUE, 
                                  WAIC = TRUE)
nimbleOptions(showCompilerOutput = FALSE)

### Summarise across MCMC samples to obtained BYM model estimates --------------

## Extract alpha, beta, sigma2.s and sigma2.u ----------------------------------

eco_ests <- eco_samples$summary$all.chains[c("alpha",
                                             "beta",
                                             "sigma2.s",
                                             "sigma2.u"),
                                           c("Median",
                                             "95%CI_low",
                                             "95%CI_upp")] %>%
  
  # Convert matrix to data frame
  
  base::as.data.frame() %>%
  
  # Rename columns
  
  stats::setNames(c("median", "lower_95_ci", "upper_95_ci")) %>%
  
  # Make the data frame row names an explicit variable called *parameter*
  
  tibble::rownames_to_column(var = "parameter") %>%
  
  # Exponentiate the parameter and 95% CI estimates for alpha to obtain overall SIR effect measure
  
  dplyr::mutate(dplyr::across(c(median, lower_95_ci, upper_95_ci), ~ base::ifelse(parameter %in% c("alpha", "beta"), base::exp(.), NA), .names = "{.col}_exp")) %>%
  
  # Format estimates as text variables with three decimal places
  
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ base::ifelse(base::is.na(.), "", base::format(base::round(., 3), nsmall = 3))))

### Print model estimates and WIAC ---------------------------------------------

base::print(list("Model estimate" = eco_ests,
                 "WAIC" = eco_samples$WAIC$WAIC))

### Summarise across MCMC samples to obtained smoothed SIR estimates -----------

# Create data set of smoothed and residual risk of lip cancer

eco_sf <- lip_sf %>%
  
  # Smoothed SIR estimate
  
  dplyr::bind_cols(eco_samples$summary$all.chains[base::paste0("smoothed.sir[", 1:N, "]"), c("Median", "95%CI_low", "95%CI_upp")] %>%
                     base::as.data.frame() %>%
                     stats::setNames(c("eco_smoothed_sir_est", "eco_smoothed_sir_l95", "eco_smoothed_sir_uci"))) %>%
  
  # Residual SIR estimate
  
  dplyr::bind_cols(eco_samples$summary$all.chains[base::paste0("residual.sir[", 1:N, "]"), c("Median", "95%CI_low", "95%CI_upp")] %>%
                     base::as.data.frame() %>%
                     stats::setNames(c("eco_residual_sir_est", "eco_residual_sir_l95", "eco_residual_sir_uci"))) %>%
  
  # Residual s SIR estimate
  
  dplyr::bind_cols(eco_samples$summary$all.chains[base::paste0("residual.sir.s[", 1:N, "]"), c("Median", "95%CI_low", "95%CI_upp")] %>%
                     base::as.data.frame() %>%
                     stats::setNames(c("eco_residual_sir_s_est", "eco_residual_sir_s_l95", "eco_residual_sir_s_uci"))) %>%
  
  # Residual u SIR estimate
  
  dplyr::bind_cols(eco_samples$summary$all.chains[base::paste0("residual.sir.u[", 1:N, "]"), c("Median", "95%CI_low", "95%CI_upp")] %>%
                     base::as.data.frame() %>%
                     stats::setNames(c("eco_residual_sir_u_est", "eco_residual_sir_u_l95", "eco_residual_sir_u_uci"))) %>%
  
  # Posterior probabilities for residual risk > 1
  
  dplyr::bind_cols(eco_samples$summary$all.chains[base::paste0("residual.sir.gt1[", 1:N, "]"), "Mean"] %>%
                     base::as.data.frame() %>%
                     stats::setNames("eco_residual_sir_gt1")) %>%
  
  # Posterior probabilities for residual risk < 1
  
  dplyr::bind_cols(eco_samples$summary$all.chains[base::paste0("residual.sir.lt1[", 1:N, "]"), "Mean"] %>%
                     base::as.data.frame() %>%
                     stats::setNames("eco_residual_sir_lt1")) %>%
  
  # Annotate exceedence probabilities
  
  dplyr::mutate(eco_residual_sir_cat = dplyr::case_when(eco_residual_sir_gt1 > 0.995 ~ "++",
                                                        eco_residual_sir_gt1 > 0.975 ~ "+",
                                                        eco_residual_sir_lt1 > 0.995 ~ "--",
                                                        eco_residual_sir_lt1 > 0.975 ~ "-",
                                                        .default = ""))

### Map smoothed SIR adjusted for percent population in AFF --------------------

map_eco_residual_sir <- ggplot() +
  
  # Add SIR estimates
  
  geom_sf(data = eco_sf,
          aes(fill = eco_residual_sir_est),
          col = "grey50",
          linewidth = 0.1) +
  
  # Add SIR posterior probabilities
  
  ggplot2::geom_point(data = eco_sf %>% dplyr::filter(eco_residual_sir_cat != ""),
                      mapping = ggplot2::aes(shape = eco_residual_sir_cat,
                                             geometry = geom),
                      fill = "black",
                      stat = "sf_coordinates") +
  
  # Add Scotland boundary
  
  ggplot2::geom_sf(data = bnd_sf,
                   fill = NA,
                   colour = "Black",
                   linewidth = 0.25) +
  
  # Add fill aesthetic gradient mapping
  
  ggplot2::scale_fill_gradient2("AFF adjusted SIR",
                                low = col_pal["low"], # blue
                                mid = col_pal["mid"], # yellow
                                high = col_pal["high"], # darkblue
                                midpoint = 0,
                                na.value = NA,
                                limits = c(base::exp(-max_log_sir), base::exp(max_log_sir)),
                                breaks = base::exp(base::seq(-max_log_sir, max_log_sir, 0.5)),
                                labels = base::format(base::round(exp(seq(-max_log_sir, max_log_sir, 0.5)), 2), nsmall = 2),
                                trans = "log") +
  
  # Add shape aesthic mapping
  
  ggplot2::scale_shape_manual("Posterior probabilities",
                              values = c("++" = 24,
                                         "+" = 2,
                                         "-" = 6,
                                         "--" = 25),
                              breaks = c("++", "+", "-", "--"),
                              labels = c("++" = "SIR > 1, p > 97.5%",
                                         "+" = "SIR > 1, p > 95.0%",
                                         "-" = "SIR < 1, p > 95.0%",
                                         "--" = "SIR < 1, p > 97.5%")) +
  
  # Format legends
  
  guides(fill = guide_colourbar(frame.colour = "black", 
                                frame.linewidth = 0.25,
                                ticks.colour = "black",
                                ticks.linewidth = 0.25,
                                order = 1),
         shape = guide_legend(order = 2)) +
  
  # Format map
  
  theme_bw() + 
  
  theme(axis.title = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        
        legend.position = "inside",
        legend.margin = ggplot2::margin(0, 0, 0, 0),
        legend.position.inside = c(0.15, 0.88),
        
        panel.grid = ggplot2::element_blank(),
        
        text = ggplot2::element_text(size = 12))

## Output observed, BYM smoothed and BYM adjusted SIR --------------------------

map_obs_bym_eco <- (map_obs_sir | map_bym_residual_sir | map_eco_residual_sir) +
  patchwork::plot_annotation(title = "Observed, BYM smoothed and AFF adjusted SIR for male lip cancer in Scotland (1975-1980)", 
                             theme = theme(plot.title = element_text(hjust = 0.5)))

base::print(map_obs_bym_eco)

# Save map - not run
# 
# ggplot2::ggsave(map_obs_bym_eco,
#                 filename = "Figures/Observed, BYM smoothed and AFF adjusted SIR for male lip cancer in Scotland (1975-1980).png",
#                 device = "png",
#                 height = 99.3 * 4,
#                 width = 228.6 * 4,
#                 units = "mm",
#                 bg = "white",
#                 dpi = 300)