#!/usr/bin/env Rscript

# This script tests the following:
# Pagel's lambda
# Phylo only
# Spatial only
# Dual model

suppressPackageStartupMessages({
  library(ape)
  library(geiger)
  library(INLA)
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

### Data ####
data = read.csv(filename)

tree = read.tree("jaeger_pruned.tree")

#### Pagel's Lambda ####
y = factor(data$y)
names(y) = data$Language_ID
pagels_lambda = fitDiscrete(tree, 
                            y, 
                            transform = "lambda")

#### INLA Set up ####
kappa = 2 # smoothness parameter as recommended by Dinnage et al. (2020)
sigma = c(1, 1.15) # Sigma parameter. First value is not used. 

# Phylogenetic Precision matrix
phylo_covar_mat <- ape::vcv(tree)
phylo_covar_mat <- phylo_covar_mat / max(phylo_covar_mat)
# The diagonal of phylo_covar_mat should inform our prior
phylo_prec_mat = cov2precision(phylo_covar_mat)

# Spatial Precision matrix
spatial_covar_mat = varcov.spatial(data[,c("Longitude", "Latitude")], 
                                   cov.pars = sigma, 
                                   kappa = kappa)$varcov
dimnames(spatial_covar_mat) = list(
  data$Language_ID, data$Language_ID
)
spatial_prec_mat = cov2precision(spatial_covar_mat)


pcprior = list(prec = list(
  prior="pc.prec",
  param = c(1, 0.1)) # probability that lambda is 0.1 is 10%
)

data$Language_ID2 = data$Language_ID
data$error_id = data$Language_ID
rownames(data) = data$Language_ID

#### Phylo Only model ####
lambdaonly_model = inla(formula = y ~
                      f(Language_ID,
                        model = "generic0",
                        Cmatrix = phylo_prec_mat,
                        constr = TRUE,
                        hyper = pcprior,
                        values = rownames(phylo_prec_mat)) +
                      f(Language_ID2,
                        model = "iid",
                        hyper = pcprior,
                        constr = TRUE),
                    family = "binomial",
                    data = data)

#### Spatial Only model ####
spatialonly_model = inla(formula = y ~
                           f(Language_ID,
                             model = "generic0",
                             Cmatrix = spatial_prec_mat,
                             constr = TRUE,
                             hyper = pcprior,
                             values = rownames(spatial_prec_mat)) +
                           f(Language_ID2,
                             model = "iid",
                             hyper = pcprior,
                             constr = TRUE),
                         family = "binomial",
                         data = data)

#### Dual Model ####
dual_model = inla(y ~ 
                    f(Language_ID,
                      model = "generic0",
                      Cmatrix = phylo_prec_mat,
                      constr = TRUE,
                      hyper = pcprior,
                      values = rownames(phylo_prec_mat)) + 
                    f(Language_ID2,
                      model = "generic0",
                      Cmatrix = spatial_prec_mat,
                      constr = TRUE,
                      hyper = pcprior,
                      values = rownames(spatial_prec_mat)) +
                    f(error_id,
                      model = "iid",
                      constr = TRUE,
                      hyper = pcprior) , 
                  data = data)


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