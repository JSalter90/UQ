# Experiment using each of the left-out runs as the 'truth' in turn
# and seeing how well we can identify this

# Define observations using each of the 250 left out runs
# Calculate implausibility using different % error
# Assess how often the different metrics rule out the truth
# This might not be helpful if NROY = 100%
# So also calculate across other locations, and find ranking of this point - does it minimise impl?
val_inds <- readRDS("papers/raikoke/data/val_inds.rds")
tDataT3 <- readRDS("papers/raikoke/data/tDataT3.rds")

#### Also need predictions over the ensemble somewhere ####

tmp_data <- tDataT3[val_inds,] # selecting the 250 left-out runs
tmp_inds <- val_inds


pseudo_var <- k * obs_var[1]

PseudoExperiment <- function(val_data, val_inds, em_pred, obs_error, ){
  
  overall_impl <- numeric(250) # implausibility of assumed obs in each experiment, for overall emulator
  total_matches <- numeric(250) # how many METs are considered not implausible, in each experiment
  overall_size <- pseudo_size <- cons_size <- numeric(250) # size of the different spaces
  
  n <- nrow(val_data)
  
  for (i in 1:n){
    pseudo_obs <- val_data$LogTotal[i]

    # Calculate overall implausibility
    overall_impl[i] <- abs(em_pred$overall$Mean[val_inds[i]] - pseudo_obs) / sqrt(em_pred$overall$SD[val_inds[i]]^2 + obs_error)

    # Size of overall NROY space
    impl_all <- (abs(em_pred$overall$Mean - pseudo_obs) / sqrt(em_pred$overall$SD^2 + obs_error))[val_inds]
    overall_size[i] <- sum(impl_all < 3) / n 
    
    # Calculate implausibility for each MET
    impl_MET <- matrix(0, n, 18)
    for (j in 1:18){
      impl_MET[,j] <- abs(em_pred$met[[j]]$Mean[val_inds] - pseudo_obs) / sqrt(em_pred$met[[j]]$SD[val_inds]^2 + obs_error)
    }
    
    total_matches[i] <- sum(impl_MET[i,] < bound) # count how many METs the chosen 'obs' are in NROY for 

    # Size of pseudo NROY space
    pseudo_size[i] <- sum(apply(impl_MET < 3, 1, sum) > 8) / n
    
    # Size of conservative NROY space
    cons_size[i] <- sum(apply(impl_MET < 3, 1, sum) > 0) / n
  }
  
  return(list())
}

# How often rule out the truth using overall emulator:
sum(pseudo_impl > 3)

# How often rule out the truth using conservative definition:
sum(pseudo_impl_matches == 0)

# How often rule out the truth using random choice of MET:
sum(pseudo_impl_matches < 9)

# Sizes of the different spaces, across the 250 experiments:
summary(overall_size)
summary(pseudo_size)
summary(cons_size)
