# Define design for evaluating emulators with
N <- 10^5
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

# Predict at T3/T5/T7 using a) overall emulator and b) each of the 18 MET emulators
ObsPred_T3 <- PredictBoth(Em_Obs_T3, Em_Obs_T3_MET, CalDesignScaled)
ObsPred_T5 <- PredictBoth(Em_Obs_T5, Em_Obs_T5_MET, CalDesignScaled)
ObsPred_T7 <- PredictBoth(Em_Obs_T7, Em_Obs_T7_MET, CalDesignScaled)

# Calculating implausibility for the overall emulator
impl_T3 <- abs(scale_output + ObsPred_T3$overall$Mean - obs[1]) / sqrt(ObsPred_T3$overall$SD^2 + obs_var[1])
impl_T5 <- abs(scale_output + ObsPred_T5$overall$Mean - obs[2]) / sqrt(ObsPred_T5$overall$SD^2 + obs_var[2])
impl_T7 <- abs(scale_output + ObsPred_T7$overall$Mean - obs[3]) / sqrt(ObsPred_T7$overall$SD^2 + obs_var[3])

bound <- 3 # standard threshold
c(sum(impl_T3 < bound), sum(impl_T5 < bound), sum(impl_T7 < bound)) # NROY at each time point independently
sum(impl_T3 < bound & impl_T5 < bound & impl_T7 < bound) # NROY for all time points

# Repeating for each MET emulator, each time
impl_T3_MET <- matrix(0, N, 18)
for (j in 1:18){
  impl_T3_MET[,j] <- abs(scale_output + ObsPred_T3$met[[j]]$Mean - obs[1]) / sqrt(ObsPred_T3$met[[j]]$SD^2 + obs_var[1])
}

impl_T5_MET <- matrix(0, N, 18)
for (j in 1:18){
  impl_T5_MET[,j] <- abs(scale_output + ObsPred_T5$met[[j]]$Mean - obs[2]) / sqrt(ObsPred_T5$met[[j]]$SD^2 + obs_var[2])
}

impl_T7_MET <- matrix(0, N, 18)
for (j in 1:18){
  impl_T7_MET[,j] <- abs(scale_output + ObsPred_T7$met[[j]]$Mean - obs[3]) / sqrt(ObsPred_T7$met[[j]]$SD^2 + obs_var[3])
}

# Finding size of NROY under different assumptions, for different time points, different emulators
NROY_T3 <- data.frame(Overall = impl_T3 < bound, # just considering the overall emulator
                      Strict = apply(impl_T3_MET < bound, 1, sum) == 18, # all 18 METs must pass
                      Conservative = apply(impl_T3_MET < bound, 1, sum) > 0, # any single MET can be not implausible
                      Pseudo = apply(impl_T3_MET < bound, 1, sum) >= 9) # at least half of METs are not implausible

# Same, but for T3 and T5
NROY_T5 <- data.frame(Overall = impl_T3 < bound & impl_T5 < bound,
                      Strict = apply(impl_T3_MET < bound, 1, sum) == 18 & apply(impl_T5_MET < bound, 1, sum) == 18,
                      Conservative = apply(impl_T3_MET < bound, 1, sum) > 0 & apply(impl_T5_MET < bound, 1, sum) > 0,
                      Pseudo = apply(impl_T3_MET < bound, 1, sum) >= 9 & apply(impl_T5_MET < bound, 1, sum) >= 9)

# Same, but for T3, T5 and T7
NROY_T7 <- data.frame(Overall = impl_T3 < bound & impl_T5 < bound & impl_T7 < bound,
                      Strict = apply(impl_T3_MET < bound, 1, sum) == 18 & apply(impl_T5_MET < bound, 1, sum) == 18 & apply(impl_T7_MET < bound, 1, sum) == 18,
                      Conservative = apply(impl_T3_MET < bound, 1, sum) > 0 & apply(impl_T5_MET < bound, 1, sum) > 0 & apply(impl_T7_MET < bound, 1, sum) > 0,
                      Pseudo = apply(impl_T3_MET < bound, 1, sum) >= 9 & apply(impl_T5_MET < bound, 1, sum) >= 9 & apply(impl_T7_MET < bound, 1, sum) >= 9)



