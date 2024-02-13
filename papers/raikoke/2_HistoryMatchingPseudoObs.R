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

Exp_T3_var1 <- PseudoExperiment(tDataT3, val_inds, EnsPred_T3, obs_error = 1 * obs_var[1])
Exp_T3_var01 <- PseudoExperiment(tDataT3, val_inds, EnsPred_T3, obs_error = 0.1 * obs_var[1])
Exp_T3_var001 <- PseudoExperiment(tDataT3, val_inds, EnsPred_T3, obs_error = 0.01 * obs_var[1])

# How often rule out the truth using overall emulator:
sum(Exp_T3_var1$overall_impl > 3)

# How often rule out the truth using conservative definition:
sum(Exp_T3_var1$total_matches == 0)

# How often rule out the truth using random choice of MET:
sum(Exp_T3_var1$total_matches < 9)

# Sizes of the different spaces, across the 250 experiments:
summary(Exp_T3_var1$overall_size)
summary(Exp_T3_var1$pseudo_size)
summary(Exp_T3_var1$cons_size)

# Or combine
data.frame(Type = c('Pseudo', 'Overall', 'Cons'),
           Errors = c(sum(Exp_T3_var1$total_matches < 9),
                      sum(Exp_T3_var1$overall_impl > 3), 
                      sum(Exp_T3_var1$total_matches == 0)))






