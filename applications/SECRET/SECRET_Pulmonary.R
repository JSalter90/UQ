# SECRET competition, September 2023
# Code for pulmonary model, Exeter team

# See Example_EmulateTS.html for more detail on the emulation approach, applied to a single time series (generalises to other high dimensional output, e.g., stacking together multiple series as here)

library(R.matlab)

# Basis emulation functions are found at https://github.com/JSalter90/UQ
# Also reads a file from https://github.com/BayesExeter/ExeterUQ, edit paths in Gasp.R to reflect location of this
setwd('~/Dropbox/UQ/') # edit directory on local machine to where https://github.com/JSalter90/UQ is cloned
source('code/Gasp.R')

# Load inputs
setwd('~/Dropbox/UQ/applications/SECRET')
designW1 <- readRDS('data/pulmonary/designW1.rds') # on original scale

# Scale to [-1,1]^4
designW1_em <- designW1
designW1_em$kMV <- (designW1$kMV - 9*10^4) / ((3*10^5 - 9*10^4)/2) - 1
designW1_em$alpha <- (designW1$alpha - 0.83) / ((0.89 - 0.83)/2) - 1
designW1_em$lrrA <- (designW1$lrrA - 20) / ((50 - 20)/2) - 1
designW1_em$lrrV <- (designW1$lrrV - 20) / ((50 - 20)/2) - 1

# Load outputs
outputW1 <- readRDS('data/pulmonary/outputW1.rds') # 512x3x100, ordered flow17, flow19, pressure

# Load observations
tmp_obs <- readMat('data/pulmonary/obs.mat')
obs <- cbind(tmp_obs[[1]], tmp_obs[[2]]) # 512x3, ordered flow17, flow19, pressure
rm(tmp_obs)

# At the time, the true inputs were unknown
# Included here so can plot truth alongside history matching results
truth <- data.frame(kMV = 2.5*10^5,
                    alpha = 0.885,
                    lrrA = 35,
                    lrrV = 25,
                    l = 0.1,
                    var = 5.8)


#### Wave 1 ####
# Emulate using initial 100 member LHC design
# Emulate flow17 (hereafter: flow1), flow19 (= flow2) and pressure series separately
# Should add a bit more flexibility to the emulator for each, may be helpful in identifying the error process (we need to emulate as accurately as possible everywhere in order to do so)

# Calculate basis, project
DataBasis_flow1 <- MakeDataBasis(outputW1[,1,])
q1 <- ExplainT(DataBasis_flow1, vtot = 0.95)
Coeffs_flow1 <- CalcScores(data = DataBasis_flow1$CentredField, basis = DataBasis_flow1$tBasis[,1:q1])
tData_flow1 <- data.frame(designW1_em[,1:4], Noise = runif(100), Coeffs_flow1)

DataBasis_flow2 <- MakeDataBasis(outputW1[,2,])
q2 <- ExplainT(DataBasis_flow2, vtot = 0.95)
Coeffs_flow2 <- CalcScores(data = DataBasis_flow2$CentredField, basis = DataBasis_flow2$tBasis[,1:q2])
tData_flow2 <- data.frame(designW1_em[,1:4], Noise = runif(100), Coeffs_flow2)

DataBasis_pressure <- MakeDataBasis(outputW1[,3,])
q3 <- ExplainT(DataBasis_pressure, vtot = 0.999) # different choice, as only need 3 vectors to explain 99.9% here
Coeffs_pressure <- CalcScores(data = DataBasis_pressure$CentredField, basis = DataBasis_pressure$tBasis[,1:q3])
tData_pressure <- data.frame(designW1_em[,1:4], Noise = runif(100), Coeffs_pressure)

# Emulate
em_flow1 <- BasisEmulators(tData_flow1, q1, mean_fn = 'step', maxdf = NULL, training_prop = 1)
par(mfrow=c(2,2), mar=c(4,4,2,2))
LeaveOneOut(em_flow1[[1]]);LeaveOneOut(em_flow1[[2]]);LeaveOneOut(em_flow1[[3]]);LeaveOneOut(em_flow1[[4]])

em_flow2 <- BasisEmulators(tData_flow2, q2, mean_fn = 'step', maxdf = NULL, training_prop = 1)
par(mfrow=c(2,2), mar=c(4,4,2,2))
LeaveOneOut(em_flow2[[1]]);LeaveOneOut(em_flow2[[2]]);LeaveOneOut(em_flow2[[3]]);LeaveOneOut(em_flow2[[4]])

em_pressure <- BasisEmulators(tData_pressure, q3, mean_fn = 'step', maxdf = NULL, training_prop = 1)
par(mfrow=c(2,2), mar=c(4,4,2,2))
LeaveOneOut(em_pressure[[1]]);LeaveOneOut(em_pressure[[2]]);LeaveOneOut(em_pressure[[3]])

# Predict across large LHC
BigDesign <- 2*as.data.frame(randomLHS(100000, 4)) - 1
colnames(BigDesign) <- colnames(designW1_em)[1:4]
Big_preds_flow1 <- BasisPredGasp(BigDesign, em_flow1)
Big_preds_flow2 <- BasisPredGasp(BigDesign, em_flow2)
Big_preds_pressure <- BasisPredGasp(BigDesign, em_pressure)




#### Wave 2 ####
designW2 <- readRDS('data/pulmonary/designW2.rds') # on original scale

# Scale to [-1,1]^4
designW2_em <- designW2
designW2_em$kMV <- (designW2$kMV - 9*10^4) / ((3*10^5 - 9*10^4)/2) - 1
designW2_em$alpha <- (designW2$alpha - 0.83) / ((0.89 - 0.83)/2) - 1
designW2_em$lrrA <- (designW2$lrrA - 20) / ((50 - 20)/2) - 1
designW2_em$lrrV <- (designW2$lrrV - 20) / ((50 - 20)/2) - 1

# Load outputs
outputW2 <- readRDS('data/pulmonary/outputW2.rds')



#### Wave 3 ####
designW3 <- readRDS('data/pulmonary/designW3.rds') # on original scale

# Scale to [-1,1]^4
designW3_em <- designW3
designW3_em$kMV <- (designW3$kMV - 9*10^4) / ((3*10^5 - 9*10^4)/2) - 1
designW3_em$alpha <- (designW3$alpha - 0.83) / ((0.89 - 0.83)/2) - 1
designW3_em$lrrA <- (designW3$lrrA - 20) / ((50 - 20)/2) - 1
designW3_em$lrrV <- (designW3$lrrV - 20) / ((50 - 20)/2) - 1

# Load outputs
outputW3 <- readRDS('data/pulmonary/outputW3.rds')




