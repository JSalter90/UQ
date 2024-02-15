# Load required packages, functions
source("papers/raikoke/0_Source.R")

# Define design for evaluating emulators with
N <- 10^5
# Have to generate this initially on original scale rather than [-1,1], as need to calculate MER from H, MER.F, DFAF - not a completely random design
CalDesign <- data.frame(H = runif(N, 9, 17),
                        DFAF = runif(N, 0.5, 20),
                        MER.F = runif(N, 0.33, 3),
                        rho = runif(N, 1350, 2500),
                        Length = sample(9:15, N, replace = TRUE),
                        sigU = exp(runif(N, log(0.0025), log(2.75))),
                        tauU = runif(N, 100, 900),
                        m_sigU = runif(N, 0.27, 1.74))

# MER is a function of other inputs, calculate its log
CalDesign$logMER <- log(50.7 * 10^7 * (CalDesign$H - 0.551)^(1/0.241) * CalDesign$MER.F * CalDesign$DFAF/100)

# Scale to [-1,1] for prediction
CalDesignScaled <- CalDesign
for (i in 1:ncol(CalDesignScaled)){
  ind <- which(parRanges$parameter == colnames(CalDesignScaled)[i])
  CalDesignScaled[,i] <- CalDesignScaled[,i] - parRanges$lower[ind]
  CalDesignScaled[,i] <- CalDesignScaled[,i] / ((parRanges$upper[ind] - parRanges$lower[ind])/2)
  CalDesignScaled[,i] <- CalDesignScaled[,i] - 1
}

# Load in emulators
# For the totals at T3, T5 and T7 these are available in data folder
# For others, need to fit and store locally
EmT3 <- load_ExUQmogp('papers/raikoke/data/EmT3')
EmT5 <- load_ExUQmogp('papers/raikoke/data/EmT5')
EmT7 <- load_ExUQmogp('papers/raikoke/data/EmT7')

EmT3_MET <- LoadMulti('papers/raikoke/data/EmT3_MET')
EmT5_MET <- LoadMulti('papers/raikoke/data/EmT5_MET')
EmT7_MET <- LoadMulti('papers/raikoke/data/EmT7_MET')

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

# Predict at T3/T5/T7 using a) overall emulator and b) each of the 18 MET emulators
# Each of these objects is relatively large (mean, variance at N points for 19 emulators), hence not stored in Github
PredT3 <- PredictBoth(EmT3, EmT3_MET, CalDesignScaled)
PredT5 <- PredictBoth(EmT5, EmT5_MET, CalDesignScaled)
PredT7 <- PredictBoth(EmT7, EmT7_MET, CalDesignScaled)
PredR1 <- PredictBoth(EmRegion1, EmRegion1_MET, CalDesignScaled)
PredR2 <- PredictBoth(EmRegion2, EmRegion2_MET, CalDesignScaled)
PredR3 <- PredictBoth(EmRegion3, EmRegion3_MET, CalDesignScaled)
PredR4 <- PredictBoth(EmRegion4, EmRegion4_MET, CalDesignScaled)
PredR5 <- PredictBoth(EmRegion5, EmRegion5_MET, CalDesignScaled)
PredR6 <- PredictBoth(EmRegion6, EmRegion6_MET, CalDesignScaled)
PredR7 <- PredictBoth(EmRegion7, EmRegion7_MET, CalDesignScaled)
PredR8 <- PredictBoth(EmRegion8, EmRegion8_MET, CalDesignScaled)

# Calculating implausibility for the overall emulator
implT3 <- abs(scale_output + PredT3$overall$Mean - obs[1]) / sqrt(PredT3$overall$SD^2 + obs_var[1])
implT5 <- abs(scale_output + PredT5$overall$Mean - obs[2]) / sqrt(PredT5$overall$SD^2 + obs_var[2])
implT7 <- abs(scale_output + PredT7$overall$Mean - obs[3]) / sqrt(PredT7$overall$SD^2 + obs_var[3])

bound <- 3 # standard threshold
c(sum(implT3 < bound), sum(implT5 < bound), sum(implT7 < bound)) # NROY at each time point independently
sum(implT3 < bound & implT5 < bound & implT7 < bound) # NROY for all time points

# Repeating for each MET emulator, each time
implT3_MET <- matrix(0, N, 18)
for (j in 1:18){
  implT3_MET[,j] <- abs(scale_output + PredT3$met[[j]]$Mean - obs[1]) / sqrt(PredT3$met[[j]]$SD^2 + obs_var[1])
}

implT5_MET <- matrix(0, N, 18)
for (j in 1:18){
  implT5_MET[,j] <- abs(scale_output + PredT5$met[[j]]$Mean - obs[2]) / sqrt(PredT5$met[[j]]$SD^2 + obs_var[2])
}

implT7_MET <- matrix(0, N, 18)
for (j in 1:18){
  implT7_MET[,j] <- abs(scale_output + PredT7$met[[j]]$Mean - obs[3]) / sqrt(PredT7$met[[j]]$SD^2 + obs_var[3])
}

# Finding size of NROY under different assumptions, for different time points, different emulators
NROY_T3 <- data.frame(Overall = implT3 < bound, # just considering the overall emulator
                      Strict = apply(implT3_MET < bound, 1, sum) == 18, # all 18 METs must pass
                      Conservative = apply(implT3_MET < bound, 1, sum) > 0, # any single MET can be not implausible
                      Pseudo = apply(implT3_MET < bound, 1, sum) >= 9) # at least half of METs are not implausible

# Same, but for T3 and T5
NROY_T5 <- data.frame(Overall = implT3 < bound & implT5 < bound,
                      Strict = apply(implT3_MET < bound, 1, sum) == 18 & apply(implT5_MET < bound, 1, sum) == 18,
                      Conservative = apply(implT3_MET < bound, 1, sum) > 0 & apply(implT5_MET < bound, 1, sum) > 0,
                      Pseudo = apply(implT3_MET < bound, 1, sum) >= 9 & apply(implT5_MET < bound, 1, sum) >= 9)

# Same, but for T3, T5 and T7
NROY_T7 <- data.frame(Overall = implT3 < bound & implT5 < bound & implT7 < bound,
                      Strict = apply(implT3_MET < bound, 1, sum) == 18 & apply(implT5_MET < bound, 1, sum) == 18 & apply(implT7_MET < bound, 1, sum) == 18,
                      Conservative = apply(implT3_MET < bound, 1, sum) > 0 & apply(implT5_MET < bound, 1, sum) > 0 & apply(implT7_MET < bound, 1, sum) > 0,
                      Pseudo = apply(implT3_MET < bound, 1, sum) >= 9 & apply(implT5_MET < bound, 1, sum) >= 9 & apply(implT7_MET < bound, 1, sum) >= 9)

# N+S (R1, R2) FIX OBS
implR1 <- abs(scale_output + PredR1$overall$Mean - obs[1]) / sqrt(PredR1$overall$SD^2 + obs_var[1])
implR2 <- abs(scale_output + PredR2$overall$Mean - obs[2]) / sqrt(PredR2$overall$SD^2 + obs_var[2])

implR1_MET <- matrix(0, N, 18)
for (j in 1:18){
  implR1_MET[,j] <- abs(scale_output + PredR1$met[[j]]$Mean - obs[1]) / sqrt(PredR1$met[[j]]$SD^2 + obs_var[1])
}

implR2_MET <- matrix(0, N, 18)
for (j in 1:18){
  implR2_MET[,j] <- abs(scale_output + PredR2$met[[j]]$Mean - obs[1]) / sqrt(PredR2$met[[j]]$SD^2 + obs_var[1])
}

NROY_NS <- data.frame(Overall = implR1 < bound & implR2 < bound,
                      Strict = apply(implR1_MET < bound & implR2_MET < bound, 1, sum) == 18,
                      Conservative = apply(implR1_MET < bound & implR2_MET < bound, 1, sum) > 0,
                      Pseudo = apply(implR1_MET < bound & implR2_MET < bound, 1, sum) >= 9)

# W+E (R3, R4) FIX OBS
implR3 <- abs(scale_output + PredR3$overall$Mean - obs[1]) / sqrt(PredR3$overall$SD^2 + obs_var[1])
implR4 <- abs(scale_output + PredR4$overall$Mean - obs[2]) / sqrt(PredR4$overall$SD^2 + obs_var[2])

implR3_MET <- matrix(0, N, 18)
for (j in 1:18){
  implR3_MET[,j] <- abs(scale_output + PredR3$met[[j]]$Mean - obs[1]) / sqrt(PredR3$met[[j]]$SD^2 + obs_var[1])
}

implR4_MET <- matrix(0, N, 18)
for (j in 1:18){
  implR4_MET[,j] <- abs(scale_output + PredR4$met[[j]]$Mean - obs[1]) / sqrt(PredR4$met[[j]]$SD^2 + obs_var[1])
}

NROY_NS <- data.frame(Overall = implR3 < bound & implR4 < bound,
                      Strict = apply(implR3_MET < bound & implR4_MET < bound, 1, sum) == 18,
                      Conservative = apply(implR3_MET < bound & implR4_MET < bound, 1, sum) > 0,
                      Pseudo = apply(implR3_MET < bound & implR4_MET < bound, 1, sum) >= 9)

# 4 regions
implR5 <- abs(scale_output + PredR5$overall$Mean - obs[1]) / sqrt(PredR5$overall$SD^2 + obs_var[1])
implR6 <- abs(scale_output + PredR6$overall$Mean - obs[2]) / sqrt(PredR6$overall$SD^2 + obs_var[2])
implR7 <- abs(scale_output + PredR7$overall$Mean - obs[2]) / sqrt(PredR7$overall$SD^2 + obs_var[2])
implR8 <- abs(scale_output + PredR8$overall$Mean - obs[2]) / sqrt(PredR8$overall$SD^2 + obs_var[2])

implR5_MET <- matrix(0, N, 18)
for (j in 1:18){
  implR5_MET[,j] <- abs(scale_output + PredR5$met[[j]]$Mean - obs[1]) / sqrt(PredR5$met[[j]]$SD^2 + obs_var[1])
}

implR6_MET <- matrix(0, N, 18)
for (j in 1:18){
  implR6_MET[,j] <- abs(scale_output + PredR6$met[[j]]$Mean - obs[1]) / sqrt(PredR6$met[[j]]$SD^2 + obs_var[1])
}

implR7_MET <- matrix(0, N, 18)
for (j in 1:18){
  implR7_MET[,j] <- abs(scale_output + PredR7$met[[j]]$Mean - obs[1]) / sqrt(PredR7$met[[j]]$SD^2 + obs_var[1])
}

implR8_MET <- matrix(0, N, 18)
for (j in 1:18){
  implR8_MET[,j] <- abs(scale_output + PredR8$met[[j]]$Mean - obs[1]) / sqrt(PredR8$met[[j]]$SD^2 + obs_var[1])
}

NROY_4R <- data.frame(Overall = implR5 < bound & implR6 < bound & implR7 < bound & implR8 < bound,
                      Strict = apply(implR5_MET < bound & implR6_MET < bound & implR7_MET < bound & implR8_MET < bound, 1, sum) == 18,
                      Conservative = apply(implR5_MET < bound & implR6_MET < bound & implR7_MET < bound & implR8_MET < bound, 1, sum) > 0,
                      Pseudo = apply(implR5_MET < bound & implR6_MET < bound & implR7_MET < bound & implR8_MET < bound, 1, sum) >= 9)

# Instead find 2nd max
impl4R_k2 <- apply(cbind(implR5, implR6, implR7, implR8), 1, kth_max, k = 2)
impl4R_k2_MET <- apply(abind(implR5, implR6, implR7, implR8), c(1,2), kth_max, k = 2)

NROY_4R_k2 <- data.frame(Overall = impl4R_k2 < bound,
                         Strict = apply(impl4R_k2_MET < bound, 1, sum) == 18,
                         Conservative = apply(impl4R_k2_MET < bound, 1, sum) > 0,
                         Pseudo = apply(impl4R_k2_MET < bound, 1, sum) >= 9)

