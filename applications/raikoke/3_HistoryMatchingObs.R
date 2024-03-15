# Load required packages, functions
source("applications/raikoke/0_Source.R")

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

# In general, scale down output by 0.95 when compare to observations
# All on log scale here
scale_output <- log(0.95)

# Load in emulators
# For the totals at T1, T3 and T5 these are available in data folder
# For others, need to fit and store locally
EmT1 <- load_ExUQmogp('applications/raikoke/data/EmT1/train')
EmT3 <- load_ExUQmogp('applications/raikoke/data/EmT3/train')
EmT5 <- load_ExUQmogp('applications/raikoke/data/EmT5/train')

EmT1_MET <- LoadMulti('applications/raikoke/data/EmT1_MET', 'train')
EmT3_MET <- LoadMulti('applications/raikoke/data/EmT3_MET', 'train')
EmT5_MET <- LoadMulti('applications/raikoke/data/EmT5_MET', 'train')

# Regions, overall
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

# Predict at T1/T3/T5 using a) overall emulator and b) each of the 18 MET emulators
# Each of these objects is relatively large (mean, variance at N points for 19 emulators), hence not stored in Github
PredT1 <- PredictBoth(EmT1, EmT1_MET, CalDesignScaled)
PredT3 <- PredictBoth(EmT3, EmT3_MET, CalDesignScaled)
PredT5 <- PredictBoth(EmT5, EmT5_MET, CalDesignScaled)
PredR1 <- PredictBoth(EmRegion1, EmRegion1_MET, CalDesignScaled)
PredR2 <- PredictBoth(EmRegion2, EmRegion2_MET, CalDesignScaled)
PredR3 <- PredictBoth(EmRegion3, EmRegion3_MET, CalDesignScaled)
PredR4 <- PredictBoth(EmRegion4, EmRegion4_MET, CalDesignScaled)
PredR5 <- PredictBoth(EmRegion5, EmRegion5_MET, CalDesignScaled)
PredR6 <- PredictBoth(EmRegion6, EmRegion6_MET, CalDesignScaled)
PredR7 <- PredictBoth(EmRegion7, EmRegion7_MET, CalDesignScaled)
PredR8 <- PredictBoth(EmRegion8, EmRegion8_MET, CalDesignScaled)

# Calculating implausibility for the overall emulator
implT1 <- abs(scale_output + PredT1$overall$Mean - subset(obs, Type == 'T1')$Mean) / sqrt(PredT1$overall$SD^2 + subset(obs, Type == 'T1')$Var)
implT3 <- abs(scale_output + PredT3$overall$Mean - subset(obs, Type == 'T3')$Mean) / sqrt(PredT3$overall$SD^2 + subset(obs, Type == 'T3')$Var)
implT5 <- abs(scale_output + PredT5$overall$Mean - subset(obs, Type == 'T5')$Mean) / sqrt(PredT5$overall$SD^2 + subset(obs, Type == 'T5')$Var)

bound <- 3 # standard threshold
c(sum(implT1 < bound), sum(implT3 < bound), sum(implT5 < bound)) # NROY at each time point independently
sum(implT1 < bound & implT3 < bound & implT5 < bound) # NROY for all time points

# Repeating for each MET emulator, each time
implT1_MET <- matrix(0, N, 18)
for (j in 1:18){
  implT1_MET[,j] <- abs(scale_output + PredT1$met[[j]]$Mean - subset(obs, Type == 'T1')$Mean) / sqrt(PredT1$met[[j]]$SD^2 + subset(obs, Type == 'T1')$Var)
}

implT3_MET <- matrix(0, N, 18)
for (j in 1:18){
  implT3_MET[,j] <- abs(scale_output + PredT3$met[[j]]$Mean - subset(obs, Type == 'T3')$Mean) / sqrt(PredT3$met[[j]]$SD^2 + subset(obs, Type == 'T3')$Var)
}

implT5_MET <- matrix(0, N, 18)
for (j in 1:18){
  implT5_MET[,j] <- abs(scale_output + PredT5$met[[j]]$Mean - subset(obs, Type == 'T5')$Mean) / sqrt(PredT5$met[[j]]$SD^2 + subset(obs, Type == 'T5')$Var)
}

# Also consider percentage of NAME ensemble that is ruled out
implT1Ens <- abs(scale_output + tDataT1$LogTotal - subset(obs, Type == 'T1')$Mean) / sqrt(subset(obs, Type == 'T1')$Var)
implT3Ens <- abs(scale_output + tDataT3$LogTotal - subset(obs, Type == 'T3')$Mean) / sqrt(subset(obs, Type == 'T3')$Var)
implT5Ens <- abs(scale_output + tDataT5$LogTotal - subset(obs, Type == 'T5')$Mean) / sqrt(subset(obs, Type == 'T5')$Var)

# Finding size of NROY under different assumptions, for different time points, different emulators
NROY_T1 <- data.frame(Overall = sum(implT1 < bound) / N, # just considering the overall emulator
                      Simulator = sum(implT1Ens < bound) / length(implT1Ens), # only using true output
                      Conservative = sum(apply(implT1_MET < bound, 1, sum) > 0) / N, # any single MET can be not implausible
                      Pseudo = sum(apply(implT1_MET < bound, 1, sum) >= 9) / N) # at least half of METs are not implausible

# Same, but for T1 and T3
NROY_T3 <- data.frame(Overall = sum(implT1 < bound & implT3 < bound) / N,
                      Simulator = sum(implT1Ens < bound & implT3Ens < bound) / length(implT1Ens),
                      Conservative = sum(apply(implT1_MET < bound, 1, sum) > 0 & apply(implT3_MET < bound, 1, sum) > 0) / N,
                      Pseudo = sum(apply(implT1_MET < bound, 1, sum) >= 9 & apply(implT3_MET < bound, 1, sum) >= 9) / N)

# Same, but for T1, T3 and T5
NROY_T5 <- data.frame(Overall = sum(implT1 < bound & implT3 < bound & implT5 < bound) / N,
                      Simulator = sum(implT1Ens < bound & implT3Ens < bound & implT5Ens < bound) / length(implT1Ens),
                      Conservative = sum(apply(implT1_MET < bound, 1, sum) > 0 & apply(implT3_MET < bound, 1, sum) > 0 & apply(implT5_MET < bound, 1, sum) > 0) / N,
                      Pseudo = sum(apply(implT1_MET < bound, 1, sum) >= 9 & apply(implT3_MET < bound, 1, sum) >= 9 & apply(implT5_MET < bound, 1, sum) >= 9) / N)

# N+S (R1, R2)
implR1 <- abs(scale_output + PredR1$overall$Mean - subset(obs, Type == 'N')$Mean) / sqrt(PredR1$overall$SD^2 + subset(obs, Type == 'N')$Var)
implR2 <- abs(scale_output + PredR2$overall$Mean - subset(obs, Type == 'S')$Mean) / sqrt(PredR2$overall$SD^2 + subset(obs, Type == 'S')$Var)

implR1_MET <- matrix(0, N, 18)
for (j in 1:18){
  implR1_MET[,j] <- abs(scale_output + PredR1$met[[j]]$Mean - subset(obs, Type == 'N')$Mean) / sqrt(PredR1$met[[j]]$SD^2 + subset(obs, Type == 'N')$Var)
}

implR2_MET <- matrix(0, N, 18)
for (j in 1:18){
  implR2_MET[,j] <- abs(scale_output + PredR2$met[[j]]$Mean - subset(obs, Type == 'S')$Mean) / sqrt(PredR2$met[[j]]$SD^2 + subset(obs, Type == 'S')$Var)
}

tData_regions <- readRDS("applications/raikoke/data/tData_regions.rds")
implR1Ens <- abs(scale_output + tData_regions[[1]]$LogTotal - subset(obs, Type == 'N')$Mean) / sqrt(subset(obs, Type == 'N')$Var)
implR2Ens <- abs(scale_output + tData_regions[[2]]$LogTotal - subset(obs, Type == 'S')$Mean) / sqrt(subset(obs, Type == 'S')$Var)

NROY_NS <- data.frame(Overall = sum(implR1 < bound & implR2 < bound) / N,
                      Simulator = sum(implR1Ens < bound & implR2Ens < bound) / length(implR1Ens),
                      Conservative = sum(apply(implR1_MET < bound & implR2_MET < bound, 1, sum) > 0) / N,
                      Pseudo = sum(apply(implR1_MET < bound & implR2_MET < bound, 1, sum) >= 9) / N)

# W+E (R3, R4)
implR3 <- abs(scale_output + PredR3$overall$Mean - subset(obs, Type == 'W')$Mean) / sqrt(PredR3$overall$SD^2 + subset(obs, Type == 'W')$Var)
implR4 <- abs(scale_output + PredR4$overall$Mean - subset(obs, Type == 'E')$Mean) / sqrt(PredR4$overall$SD^2 + subset(obs, Type == 'E')$Var)

implR3_MET <- matrix(0, N, 18)
for (j in 1:18){
  implR3_MET[,j] <- abs(scale_output + PredR3$met[[j]]$Mean - subset(obs, Type == 'W')$Mean) / sqrt(PredR3$met[[j]]$SD^2 + subset(obs, Type == 'W')$Var)
}

implR4_MET <- matrix(0, N, 18)
for (j in 1:18){
  implR4_MET[,j] <- abs(scale_output + PredR4$met[[j]]$Mean - subset(obs, Type == 'E')$Mean) / sqrt(PredR4$met[[j]]$SD^2 + subset(obs, Type == 'E')$Var)
}

implR3Ens <- abs(scale_output + tData_regions[[3]]$LogTotal - subset(obs, Type == 'W')$Mean) / sqrt(subset(obs, Type == 'W')$Var)
implR4Ens <- abs(scale_output + tData_regions[[4]]$LogTotal - subset(obs, Type == 'E')$Mean) / sqrt(subset(obs, Type == 'E')$Var)

NROY_NS <- data.frame(Overall = sum(implR3 < bound & implR4 < bound) / N,
                      Simulator = sum(implR3Ens < bound & implR4Ens < bound) / length(implR3Ens),
                      Conservative = sum(apply(implR3_MET < bound & implR4_MET < bound, 1, sum) > 0) / N,
                      Pseudo = sum(apply(implR3_MET < bound & implR4_MET < bound, 1, sum) >= 9) / N)

# 4 regions
implR5 <- abs(scale_output + PredR5$overall$Mean - subset(obs, Type == 'NW')$Mean) / sqrt(PredR5$overall$SD^2 + subset(obs, Type == 'NW')$Var)
implR6 <- abs(scale_output + PredR6$overall$Mean - subset(obs, Type == 'NE')$Mean) / sqrt(PredR6$overall$SD^2 + subset(obs, Type == 'NE')$Var)
implR7 <- abs(scale_output + PredR7$overall$Mean - subset(obs, Type == 'SE')$Mean) / sqrt(PredR7$overall$SD^2 + subset(obs, Type == 'SE')$Var)
implR8 <- abs(scale_output + PredR8$overall$Mean - subset(obs, Type == 'SW')$Mean) / sqrt(PredR8$overall$SD^2 + subset(obs, Type == 'SW')$Var)

implR5_MET <- matrix(0, N, 18)
for (j in 1:18){
  implR5_MET[,j] <- abs(scale_output + PredR5$met[[j]]$Mean - subset(obs, Type == 'NW')$Mean) / sqrt(PredR5$met[[j]]$SD^2 + subset(obs, Type == 'NW')$Var)
}

implR6_MET <- matrix(0, N, 18)
for (j in 1:18){
  implR6_MET[,j] <- abs(scale_output + PredR6$met[[j]]$Mean - subset(obs, Type == 'NE')$Mean) / sqrt(PredR6$met[[j]]$SD^2 + subset(obs, Type == 'NE')$Var)
}

implR7_MET <- matrix(0, N, 18)
for (j in 1:18){
  implR7_MET[,j] <- abs(scale_output + PredR7$met[[j]]$Mean - subset(obs, Type == 'SE')$Mean) / sqrt(PredR7$met[[j]]$SD^2 + subset(obs, Type == 'SE')$Var)
}

implR8_MET <- matrix(0, N, 18)
for (j in 1:18){
  implR8_MET[,j] <- abs(scale_output + PredR8$met[[j]]$Mean - subset(obs, Type == 'SW')$Mean) / sqrt(PredR8$met[[j]]$SD^2 + subset(obs, Type == 'SW')$Var)
}

implR5Ens <- abs(scale_output + tData_regions[[5]]$LogTotal - subset(obs, Type == 'NW')$Mean) / sqrt(subset(obs, Type == 'NW')$Var)
implR6Ens <- abs(scale_output + tData_regions[[6]]$LogTotal - subset(obs, Type == 'NE')$Mean) / sqrt(subset(obs, Type == 'NE')$Var)
implR7Ens <- abs(scale_output + tData_regions[[7]]$LogTotal - subset(obs, Type == 'SE')$Mean) / sqrt(subset(obs, Type == 'SE')$Var)
implR8Ens <- abs(scale_output + tData_regions[[8]]$LogTotal - subset(obs, Type == 'SW')$Mean) / sqrt(subset(obs, Type == 'SW')$Var)

NROY_4R <- data.frame(Overall = sum(implR5 < bound & implR6 < bound & implR7 < bound & implR8 < bound) / N,
                      Simulator = sum(implR5Ens < bound & implR6Ens < bound & implR7Ens < bound & implR8Ens < bound) / length(implR5Ens),
                      Conservative = sum(apply(implR5_MET < bound & implR6_MET < bound & implR7_MET < bound & implR8_MET < bound, 1, sum) > 0) / N,
                      Pseudo = sum(apply(implR5_MET < bound & implR6_MET < bound & implR7_MET < bound & implR8_MET < bound, 1, sum) >= 9) / N)

# Instead find 2nd max
impl4R_k2 <- apply(cbind(implR5, implR6, implR7, implR8), 1, kth_max, k = 2)
impl4R_k2_MET <- apply(abind(implR5_MET, implR6_MET, implR7_MET, implR8_MET), c(1,2), kth_max, k = 2)
impl4R_k2_Ens <- apply(cbind(implR5Ens, implR6Ens, implR7Ens, implR8Ens), 1, kth_max, k = 2)

NROY_4R_k2 <- data.frame(Overall = sum(impl4R_k2 < bound) / N,
                         Simulator = sum(impl4R_k2_Ens < bound) / length(implR5Ens),
                         Conservative = sum(apply(impl4R_k2_MET < bound, 1, sum) > 0) / N,
                         Pseudo = sum(apply(impl4R_k2_MET < bound, 1, sum) >= 9) / N)