# Experiment using each of the left-out runs as the 'truth' in turn
# and seeing how well we can identify this

# Define observations using each of the 250 left out runs
# Calculate implausibility using different % error
# Assess how often the different metrics rule out the truth
# This might not be helpful if NROY = 100%
# So also calculate across other locations, and find ranking of this point - does it minimise impl?
val_inds <- readRDS("papers/raikoke/data/val_inds.rds")
tDataT3 <- readRDS("papers/raikoke/data/tDataT3.rds")

# Load in emulator predictions across the 1000 ensemble members (only need the 250 val points)
EnsPredT3 <- readRDS("papers/raikoke/data/EnsPredT3.rds")

Exp_T3_var1 <- PseudoExperiment(tDataT3, val_inds, EnsPredT3, obs_error = 1 * obs_var[1])
Exp_T3_var01 <- PseudoExperiment(tDataT3, val_inds, EnsPredT3, obs_error = 0.1 * obs_var[1])
Exp_T3_var001 <- PseudoExperiment(tDataT3, val_inds, EnsPredT3, obs_error = 0.01 * obs_var[1])

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

#### Repeat for other time points, combinations of regions ####
# T5
tDataT5 <- readRDS("papers/raikoke/data/tDataT5.rds")
EnsPredT5 <- readRDS("papers/raikoke/data/EnsPredT5.rds")
Exp_T5_var1 <- PseudoExperiment(tDataT5, val_inds, EnsPredT5, obs_error = 1 * obs_var[2])
Exp_T5_var01 <- PseudoExperiment(tDataT5, val_inds, EnsPredT5, obs_error = 0.1 * obs_var[2])
Exp_T5_var001 <- PseudoExperiment(tDataT5, val_inds, EnsPredT5, obs_error = 0.01 * obs_var[2])

# T7
tDataT7 <- readRDS("papers/raikoke/data/tDataT7.rds")
EnsPredT7 <- readRDS("papers/raikoke/data/EnsPredT7.rds")
Exp_T7_var1 <- PseudoExperiment(tDataT7, val_inds, EnsPredT7, obs_error = 1 * obs_var[3])
Exp_T7_var01 <- PseudoExperiment(tDataT7, val_inds, EnsPredT7, obs_error = 0.1 * obs_var[3])
Exp_T7_var001 <- PseudoExperiment(tDataT7, val_inds, EnsPredT7, obs_error = 0.01 * obs_var[3])

# Regions
tData_regions <- readRDS("papers/raikoke/data/tData_regions.rds")
# Need to load in emulators, create predictions across the ensemble
# Overall
EmRegion1 <- load_ExUQmogp('papers/raikoke/data/EmRegion/region1')
EmRegion2 <- load_ExUQmogp('papers/raikoke/data/EmRegion/region2')
EmRegion3 <- load_ExUQmogp('papers/raikoke/data/EmRegion/region3')
EmRegion4 <- load_ExUQmogp('papers/raikoke/data/EmRegion/region4')
EmRegion5 <- load_ExUQmogp('papers/raikoke/data/EmRegion/region5')
EmRegion6 <- load_ExUQmogp('papers/raikoke/data/EmRegion/region6')
EmRegion7 <- load_ExUQmogp('papers/raikoke/data/EmRegion/region7')
EmRegion8 <- load_ExUQmogp('papers/raikoke/data/EmRegion/region8')

# By MET
EmRegion1_MET <- LoadMulti('papers/raikoke/data/EmRegions_MET/region1', 'region1')
EmRegion2_MET <- LoadMulti('papers/raikoke/data/EmRegions_MET/region2', 'region2')
EmRegion3_MET <- LoadMulti('papers/raikoke/data/EmRegions_MET/region3', 'region3')
EmRegion4_MET <- LoadMulti('papers/raikoke/data/EmRegions_MET/region4', 'region4')
EmRegion5_MET <- LoadMulti('papers/raikoke/data/EmRegions_MET/region5', 'region5')
EmRegion6_MET <- LoadMulti('papers/raikoke/data/EmRegions_MET/region6', 'region6')
EmRegion7_MET <- LoadMulti('papers/raikoke/data/EmRegions_MET/region7', 'region7')
EmRegion8_MET <- LoadMulti('papers/raikoke/data/EmRegions_MET/region8', 'region8')

# Predictions
EnsPredRegion1 <- PredictBoth(EmRegion1, EmRegion1_MET, ScaledDesign)
EnsPredRegion2 <- PredictBoth(EmRegion2, EmRegion2_MET, ScaledDesign)
EnsPredRegion3 <- PredictBoth(EmRegion3, EmRegion3_MET, ScaledDesign)
EnsPredRegion4 <- PredictBoth(EmRegion4, EmRegion4_MET, ScaledDesign)
EnsPredRegion5 <- PredictBoth(EmRegion5, EmRegion5_MET, ScaledDesign)
EnsPredRegion6 <- PredictBoth(EmRegion6, EmRegion6_MET, ScaledDesign)
EnsPredRegion7 <- PredictBoth(EmRegion7, EmRegion7_MET, ScaledDesign)
EnsPredRegion8 <- PredictBoth(EmRegion8, EmRegion8_MET, ScaledDesign)

# N+S (R1+R2)



# W+E (R3+R4)

# 4 regions (R5+R6+R7+R8)

# 3/4 regions

