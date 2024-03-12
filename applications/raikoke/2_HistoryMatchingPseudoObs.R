# Experiment using each of the left-out runs as the 'truth' in turn
# and seeing how well we can identify this

# Load required packages, functions
source("applications/raikoke/0_Source.R")

# Define observations using each of the 250 left out runs
# Calculate implausibility using different % error
# Assess how often the different metrics rule out the truth
# This might not be helpful if NROY = 100%
# So also calculate across other locations, and find ranking of this point - does it minimise impl?
val_inds <- readRDS("applications/raikoke/data/val_inds.rds")
tDataT1 <- readRDS("applications/raikoke/data/tDataT1.rds")

# Load in emulator predictions across the 1000 ensemble members (only need the 250 val points)
EnsPredT1 <- readRDS("applications/raikoke/data/EnsPredT1.rds")
ObsVar <- subset(obs, Type == 'T1')$Var
Exp_T1_var1 <- PseudoExperiment(tDataT1, val_inds, EnsPredT1, obs_error = 1 * ObsVar)
Exp_T1_var01 <- PseudoExperiment(tDataT1, val_inds, EnsPredT1, obs_error = 0.1 * ObsVar)
Exp_T1_var001 <- PseudoExperiment(tDataT1, val_inds, EnsPredT1, obs_error = 0.01 * ObsVar)

# How often rule out the truth using overall emulator:
sum(Exp_T1_var1$overall_impl > 3)

# How often rule out the truth using conservative definition:
sum(Exp_T1_var1$total_matches == 0)

# How often rule out the truth using random choice of MET:
sum(Exp_T1_var1$total_matches < 9)

# Sizes of the different spaces, across the 250 experiments:
summary(Exp_T1_var1$overall_size)
summary(Exp_T1_var1$pseudo_size)
summary(Exp_T1_var1$cons_size)

# Or combine
SummariseExperiment(Exp_T1_var1)
SummariseExperiment(Exp_T1_var01)
SummariseExperiment(Exp_T1_var001)

#### Repeat for other time points, combinations of regions ####
# T3
tDataT3 <- readRDS("applications/raikoke/data/tDataT3.rds")
EnsPredT3 <- readRDS("applications/raikoke/data/EnsPredT3.rds")
ObsVar <- subset(obs, Type == 'T3')$Var
Exp_T3_var1 <- PseudoExperiment(tDataT3, val_inds, EnsPredT3, obs_error = 1 * ObsVar)
Exp_T3_var01 <- PseudoExperiment(tDataT3, val_inds, EnsPredT3, obs_error = 0.1 * ObsVar)
Exp_T3_var001 <- PseudoExperiment(tDataT3, val_inds, EnsPredT3, obs_error = 0.01 * ObsVar)

SummariseExperiment(Exp_T3_var1)
SummariseExperiment(Exp_T3_var01)
SummariseExperiment(Exp_T3_var001)

# T5
tDataT5 <- readRDS("applications/raikoke/data/tDataT5.rds")
EnsPredT5 <- readRDS("applications/raikoke/data/EnsPredT5.rds")
ObsVar <- subset(obs, Type == 'T5')$Var
Exp_T5_var1 <- PseudoExperiment(tDataT5, val_inds, EnsPredT5, obs_error = 1 * ObsVar)
Exp_T5_var01 <- PseudoExperiment(tDataT5, val_inds, EnsPredT5, obs_error = 0.1 * ObsVar)
Exp_T5_var001 <- PseudoExperiment(tDataT5, val_inds, EnsPredT5, obs_error = 0.01 * ObsVar)

SummariseExperiment(Exp_T5_var1)
SummariseExperiment(Exp_T5_var01)
SummariseExperiment(Exp_T5_var001)

# Regions
tData_regions <- readRDS("applications/raikoke/data/tData_regions.rds")
# Need to load in emulators, create predictions across the ensemble
# Overall
EmRegion1 <- load_ExUQmogp('applications/raikoke/data/EmRegion/region1')
EmRegion2 <- load_ExUQmogp('applications/raikoke/data/EmRegion/region2')
EmRegion3 <- load_ExUQmogp('applications/raikoke/data/EmRegion/region3')
EmRegion4 <- load_ExUQmogp('applications/raikoke/data/EmRegion/region4')
EmRegion5 <- load_ExUQmogp('applications/raikoke/data/EmRegion/region5')
EmRegion6 <- load_ExUQmogp('applications/raikoke/data/EmRegion/region6')
EmRegion7 <- load_ExUQmogp('applications/raikoke/data/EmRegion/region7')
EmRegion8 <- load_ExUQmogp('applications/raikoke/data/EmRegion/region8')

# By MET
EmRegion1_MET <- LoadMulti('applications/raikoke/data/EmRegion_MET/region1', 'region1')
EmRegion2_MET <- LoadMulti('applications/raikoke/data/EmRegion_MET/region2', 'region2')
EmRegion3_MET <- LoadMulti('applications/raikoke/data/EmRegion_MET/region3', 'region3')
EmRegion4_MET <- LoadMulti('applications/raikoke/data/EmRegion_MET/region4', 'region4')
EmRegion5_MET <- LoadMulti('applications/raikoke/data/EmRegion_MET/region5', 'region5')
EmRegion6_MET <- LoadMulti('applications/raikoke/data/EmRegion_MET/region6', 'region6')
EmRegion7_MET <- LoadMulti('applications/raikoke/data/EmRegion_MET/region7', 'region7')
EmRegion8_MET <- LoadMulti('applications/raikoke/data/EmRegion_MET/region8', 'region8')

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
ObsVar <- c(subset(obs, Type == 'N')$Var, subset(obs, Type == 'S')$Var)
Exp_NS_var1 <- PseudoExperiment(list(tData_regions[[1]], tData_regions[[2]]),
                                list(val_inds, val_inds),
                                list(EnsPredRegion1, EnsPredRegion2),
                                obs_error = 1 * ObsVar)
Exp_NS_var01 <- PseudoExperiment(list(tData_regions[[1]], tData_regions[[2]]),
                                list(val_inds, val_inds),
                                list(EnsPredRegion1, EnsPredRegion2),
                                obs_error = 0.1 * ObsVar)
Exp_NS_var001 <- PseudoExperiment(list(tData_regions[[1]], tData_regions[[2]]),
                                list(val_inds, val_inds),
                                list(EnsPredRegion1, EnsPredRegion2),
                                obs_error = 0.01 * ObsVar)

SummariseExperiment(Exp_NS_var1)
SummariseExperiment(Exp_NS_var01)
SummariseExperiment(Exp_NS_var001)

# W+E (R3+R4)
ObsVar <- c(subset(obs, Type == 'W')$Var, subset(obs, Type == 'E')$Var)
Exp_WE_var1 <- PseudoExperiment(list(tData_regions[[3]], tData_regions[[4]]),
                                list(val_inds, val_inds),
                                list(EnsPredRegion3, EnsPredRegion4),
                                obs_error = 1 * ObsVar)
Exp_WE_var01 <- PseudoExperiment(list(tData_regions[[3]], tData_regions[[4]]),
                                 list(val_inds, val_inds),
                                 list(EnsPredRegion3, EnsPredRegion4),
                                 obs_error = 0.1 * ObsVar)
Exp_WE_var001 <- PseudoExperiment(list(tData_regions[[3]], tData_regions[[4]]),
                                  list(val_inds, val_inds),
                                  list(EnsPredRegion3, EnsPredRegion4),
                                  obs_error = 0.01 * ObsVar)

SummariseExperiment(Exp_WE_var1)
SummariseExperiment(Exp_WE_var01)
SummariseExperiment(Exp_WE_var001)

# 4 regions (R5+R6+R7+R8)
ObsVar <- c(subset(obs, Type == 'NW')$Var, subset(obs, Type == 'NE')$Var, subset(obs, Type == 'SE')$Var, subset(obs, Type == 'SW')$Var)
Exp_4R_var1 <- PseudoExperiment(list(tData_regions[[5]], tData_regions[[6]], tData_regions[[7]], tData_regions[[8]]),
                                list(val_inds, val_inds, val_inds, val_inds),
                                list(EnsPredRegion5, EnsPredRegion6, EnsPredRegion7, EnsPredRegion8),
                                obs_error = 1 * ObsVar)
Exp_4R_var01 <- PseudoExperiment(list(tData_regions[[5]], tData_regions[[6]], tData_regions[[7]], tData_regions[[8]]),
                                 list(val_inds, val_inds, val_inds, val_inds),
                                 list(EnsPredRegion5, EnsPredRegion6, EnsPredRegion7, EnsPredRegion8),
                                 obs_error = 0.1 * ObsVar)
Exp_4R_var001 <- PseudoExperiment(list(tData_regions[[5]], tData_regions[[6]], tData_regions[[7]], tData_regions[[8]]),
                                  list(val_inds, val_inds, val_inds, val_inds),
                                  list(EnsPredRegion5, EnsPredRegion6, EnsPredRegion7, EnsPredRegion8),
                                  obs_error = 0.01 * ObsVar)

SummariseExperiment(Exp_4R_var1)
SummariseExperiment(Exp_4R_var01)
SummariseExperiment(Exp_4R_var001)

# 3/4 regions, i.e. use 2nd max impl across the 4 regions (i.e. 1 allowed to 'fail')
Exp_4Rk2_var1 <- PseudoExperiment(list(tData_regions[[5]], tData_regions[[6]], tData_regions[[7]], tData_regions[[8]]),
                                list(val_inds, val_inds, val_inds, val_inds),
                                list(EnsPredRegion5, EnsPredRegion6, EnsPredRegion7, EnsPredRegion8),
                                obs_error = 1 * ObsVar,
                                kmax = 2)
Exp_4Rk2_var01 <- PseudoExperiment(list(tData_regions[[5]], tData_regions[[6]], tData_regions[[7]], tData_regions[[8]]),
                                 list(val_inds, val_inds, val_inds, val_inds),
                                 list(EnsPredRegion5, EnsPredRegion6, EnsPredRegion7, EnsPredRegion8),
                                 obs_error = 0.1 * ObsVar,
                                 kmax = 2)
Exp_4Rk2_var001 <- PseudoExperiment(list(tData_regions[[5]], tData_regions[[6]], tData_regions[[7]], tData_regions[[8]]),
                                  list(val_inds, val_inds, val_inds, val_inds),
                                  list(EnsPredRegion5, EnsPredRegion6, EnsPredRegion7, EnsPredRegion8),
                                  obs_error = 0.01 * ObsVar,
                                  kmax = 2)

SummariseExperiment(Exp_4Rk2_var1)
SummariseExperiment(Exp_4Rk2_var01)
SummariseExperiment(Exp_4Rk2_var001)
