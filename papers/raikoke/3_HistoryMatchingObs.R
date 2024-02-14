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

# Predict at T3/T5/T7 using a) overall emulator and b) each of the 18 MET emulators
ObsPredT3 <- PredictBoth(EmT3, EmT3_MET, CalDesignScaled)
ObsPredT5 <- PredictBoth(EmT5, EmT5_MET, CalDesignScaled)
ObsPredT7 <- PredictBoth(EmT7, EmT7_MET, CalDesignScaled)

# Calculating implausibility for the overall emulator
implT3 <- abs(scale_output + ObsPredT3$overall$Mean - obs[1]) / sqrt(ObsPredT3$overall$SD^2 + obs_var[1])
implT5 <- abs(scale_output + ObsPredT5$overall$Mean - obs[2]) / sqrt(ObsPredT5$overall$SD^2 + obs_var[2])
implT7 <- abs(scale_output + ObsPredT7$overall$Mean - obs[3]) / sqrt(ObsPredT7$overall$SD^2 + obs_var[3])

bound <- 3 # standard threshold
c(sum(implT3 < bound), sum(implT5 < bound), sum(implT7 < bound)) # NROY at each time point independently
sum(implT3 < bound & implT5 < bound & implT7 < bound) # NROY for all time points

# Repeating for each MET emulator, each time
implT3_MET <- matrix(0, N, 18)
for (j in 1:18){
  implT3_MET[,j] <- abs(scale_output + ObsPredT3$met[[j]]$Mean - obs[1]) / sqrt(ObsPredT3$met[[j]]$SD^2 + obs_var[1])
}

implT5_MET <- matrix(0, N, 18)
for (j in 1:18){
  implT5_MET[,j] <- abs(scale_output + ObsPredT5$met[[j]]$Mean - obs[2]) / sqrt(ObsPredT5$met[[j]]$SD^2 + obs_var[2])
}

implT7_MET <- matrix(0, N, 18)
for (j in 1:18){
  implT7_MET[,j] <- abs(scale_output + ObsPredT7$met[[j]]$Mean - obs[3]) / sqrt(ObsPredT7$met[[j]]$SD^2 + obs_var[3])
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



