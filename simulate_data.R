# Simulate data 
set.seed(7573673)

suppressPackageStartupMessages({
  library(ape)
  library(dplyr)
})

proportions = c(0.1, 0.25, 0.4)
# We want to get the proportions in 'proportions'
# But we are happy with proportions that are within the allowable variation
allowable_variation = 0.05 
lambda_transformations = seq(0, 1, by = 0.2)
lambda_transformations[1] = 0.01
iterations = 20

# save simulations here
simulated_location = "output/simulated_data/"
dir.create(simulated_location, recursive = TRUE, showWarnings = FALSE)

# jager tree
tree = read.tree(
  "jaeger_pruned.tree"
  )

# grambank locations
grambank_df = read.delim(
  file = "glottolog-cldf_wide_df.tsv", sep = "\t") %>% 
  dplyr::select(Language_ID, Latitude, Longitude)

for(iter in 1:iterations){
  cat("Iteration", iter, "of", max(iterations), "...\n")
  for(i in seq_along(proportions)){
    desired_proportion = proportions[i]
    cat("Searching for a proportion of", desired_proportion, "...")
    for(j in seq_along(lambda_transformations)){
      searching = TRUE
      desired_lambda = lambda_transformations[j]
      cat("\twith a lambda of", desired_lambda, "...\n")
      while(searching){
        q = matrix(c(round(runif(n = 2, min = 0, max = 1), 1)), 
                   nrow = 2,
                   ncol = 2) * c(-1, 1, 1, -1)
        y = 
          geiger::sim.char(
            geiger::rescale(tree,
                            desired_lambda,
                            model = "lambda"), 
            q, 
            model="discrete")[,1,]  
        
        observed_proportion = min(prop.table(table(y)))
        if(observed_proportion > desired_proportion - allowable_variation &
           observed_proportion < desired_proportion + allowable_variation)
          searching = FALSE
        
        out_df = data.frame(y = y - 1, # make data 0 & 1 
                   Language_ID = tree$tip.label)
        out_df = left_join(out_df, grambank_df, by = "Language_ID")
        
        filename = paste0(simulated_location, "Prop", desired_proportion,
                          "_Lambda",desired_lambda, "Iter", iter, ".csv")
        
        write.csv(out_df, filename, row.names = FALSE)
      }
    }
  }
}
