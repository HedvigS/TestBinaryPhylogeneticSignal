#!/usr/bin/env Rscript

# This script tests the following:
# Pagel's lambda
# Phylo only
# Spatial only
# Dual model

suppressPackageStartupMessages({
  library(ape)
  library(geiger)
  library(brms)
  source('varcov_spatial.R')
})

## functions
cov2precision = function(cov_mat){
  cov_mat = cov_mat / 
    exp(determinant(cov_mat)$modulus[1] /
          nrow(cov_mat))
  spatial_prec_mat = solve(cov_mat)
  spatial_prec_mat
}

args = commandArgs(trailingOnly = TRUE)
filename = args[1]
filename = "output/simulated_data/Prop0.1_Lambda0.01Iter1.csv"

### Data ####
data = read.csv(filename)

data$glottocodes = data$Language_ID
data$glottocodes2 = data$Language_ID
data$glottocodes3 = data$Language_ID

tree = read.tree("jaeger_pruned.tree")

#### Pagel's Lambda ####
y = factor(data$y)
names(y) = data$Language_ID
pagels_lambda = fitDiscrete(tree, 
                            y, 
                            transform = "lambda")

#### BRMS Set up ####
# Phylogenetic Precision matrix
phylo_covar_mat <- ape::vcv(tree)

# Spatial Precision matrix
kappa = 2 # smoothness parameter as recommended by Dinnage et al. (2020)
sigma = c(1, 1.15) # Sigma parameter. First value is not used. 
spatial_covar_mat = varcov.spatial(data[,c("Longitude", "Latitude")], 
                                   cov.pars = sigma, 
                                   kappa = kappa)$varcov
dimnames(spatial_covar_mat) = list(
  data$Language_ID, data$Language_ID
)

data$Language_ID2 = data$Language_ID
data$error_id = data$Language_ID
rownames(data) = data$Language_ID

#### Phylo Only model ####
lambdaonly_model = brm(
  y ~ 1 + 
    (1|gr(glottocodes2, cov = phylogeny)) + 
    (1|glottocodes3), 
  data = data, 
  family = bernoulli(), 
  data2 = list(phylogeny = phylo_covar_mat),
  prior = c(
    prior(normal(0, 50), "Intercept"),
    prior(student_t(3, 0, 20), "sd")
  ),
  chains = 1
)

#### Spatial Only model ####
spatialonly_model = brm(
  y ~ 1 + 
    (1|gr(glottocodes1, cov = spatial)) + 
    (1|glottocodes3), 
  data = data, 
  family = bernoulli(), 
  data2 = list(spatial = spatial_covar_mat),
  prior = c(
    prior(normal(0, 50), "Intercept"),
    prior(student_t(3, 0, 20), "sd")
  ),
  chains = 1
)

#### Dual Model ####
dual_model = brm(
  y ~ 1 + 
    (1|gr(glottocodes1, cov = spatial)) + 
    (1|gr(glottocodes2, cov = phylogeny)) + 
    (1|glottocodes3), 
  data = data, 
  family = bernoulli(), 
  data2 = list(spatial = spatial_covar_mat,
               phylogeny = phylo_covar_mat),
  prior = c(
    prior(normal(0, 50), "Intercept"),
    prior(student_t(3, 0, 20), "sd")
  ),
  chains = 1
)


#### Save output ####
dir.create("output/simulated_output/", showWarnings = FALSE)

savedsimulations_location = 
  "output/simulated_data/"

simulation_name = paste0(
  tools::file_path_sans_ext(basename(filename)),
  ".RDS"
)

save_location = paste0(savedsimulations_location, simulation_name)

saveRDS(
  list(
    pagels_lambda,
    lambdaonly_model,
    spatialonly_model,
    dual_model
    ), 
  file = save_location)