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

# History match
# Can't do as usual because the variance matrices aren't fixed/known
# We have an upper bound on the variance that we can use (if rule out for high variance, also do for lower variance, given same correlation structure)
# For correlation length, need to consider different structures
# Aim to keep runs that are not implausible across different correlation structures
time <- seq(from = 0, to = 0.85, length = 513)[-1]
Matern32 <- function(t1, t2, variance, length){
  r <- sqrt((t1 - t2)^2 / length^2)
  cov <- variance * (1 + sqrt(3)*r) * exp(-sqrt(3) * r)
  return(cov)
}

# Construct a set of possible error structures
# Here use every 0.05, but implausibility calculations are efficient so could go less coarse here
l_seq <- c(1e-10, seq(from = 0.05, to = 0.85, by = 0.05))
Sigma_e <- array(0, dim = c(512,512,length(l_seq)))
for (k in 1:length(l_seq)){
  tmp <- diag(512)
  for (i in 1:512){
    for (j in 1:512){
      tmp[i,j] <- Matern32(time[i], time[j], variance = 22.5, length = l_seq[k])
    }
  }
  tmp <- tmp + 10^-6 * diag(512)
  Sigma_e[,,k] <- tmp
}

# Calculate inverse for each (needed in HM code)
Inv <- NULL
for (k in 1:length(l_seq)){
  Inv[[k]] <- GetInverse(Sigma_e[,,k]) # using bespoke inversion code, as need to assign attributes to speed up HM
}

# History match for each output, each choice of covariance structure
impl_flow1 <- matrix(0, nrow(BigDesign), length(l_seq))
impl_flow2 <- matrix(0, nrow(BigDesign), length(l_seq))
impl_pressure <- matrix(0, nrow(BigDesign), length(l_seq))

for (k in 1:length(l_seq)){
  Big_impl_flow1 <- HistoryMatch(DataBasis_flow1, obs[,1] - DataBasis_flow1$EnsembleMean, Big_preds_flow1$Expectation, Big_preds_flow1$Variance, Error = Sigma_e[,,k], Disc = 0*diag(512), weightinv = Inv[[k]])
  impl_flow1[,k] <- Big_impl_flow1$impl
  
  Big_impl_flow2 <- HistoryMatch(DataBasis_flow2, obs[,2] - DataBasis_flow2$EnsembleMean, Big_preds_flow2$Expectation, Big_preds_flow2$Variance, Error = Sigma_e[,,k], Disc = 0*diag(512), weightinv = Inv[[k]])
  impl_flow2[,k] <- Big_impl_flow2$impl
  
  Big_impl_pressure <- HistoryMatch(DataBasis_pressure, obs[,3] - DataBasis_pressure$EnsembleMean, Big_preds_pressure$Expectation, Big_preds_pressure$Variance, Error = Sigma_e[,,k], Disc = 0*diag(512), weightinv = Inv[[k]])
  impl_pressure[,k] <- Big_impl_pressure$impl
}

apply(impl_flow1 < qchisq(0.995, 512), 2, sum)
apply(impl_flow2 < qchisq(0.995, 512), 2, sum)
apply(impl_pressure < qchisq(0.995, 512), 2, sum)


#### Add plot that shows how much rule out by choice of correlation ####




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

# Construct new basis, emulators
# Add in all of wave 1 - didn't definitively rule out any x, so still trying to emulate everything accurately
DataBasis_flow1_w2 <- MakeDataBasis(cbind(outputW1[,1,], outputW2[,1,]))
q1_w2 <- ExplainT(DataBasis_flow1_w2, vtot = 0.99)
Coeffs_flow1_w2 <- CalcScores(data = DataBasis_flow1_w2$CentredField, basis = DataBasis_flow1_w2$tBasis[,1:q1_w2])
tData_flow1_w2 <- data.frame(rbind(designW1_em[,1:4], designW2_em[,1:4]), Noise = runif(200), Coeffs_flow1_w2)

DataBasis_flow2_w2 <- MakeDataBasis(cbind(outputW1[,2,], outputW2[,2,]))
q2_w2 <- ExplainT(DataBasis_flow2_w2, vtot = 0.99)
Coeffs_flow2_w2 <- CalcScores(data = DataBasis_flow2_w2$CentredField, basis = DataBasis_flow2_w2$tBasis[,1:q2_w2])
tData_flow2_w2 <- data.frame(rbind(designW1_em[,1:4], designW2_em[,1:4]), Noise = runif(200), Coeffs_flow2_w2)

DataBasis_pressure_w2 <- MakeDataBasis(cbind(outputW1[,3,], outputW2[,3,]))
q3_w2 <- ExplainT(DataBasis_pressure_w2, vtot = 0.999)
Coeffs_pressure_w2 <- CalcScores(data = DataBasis_pressure_w2$CentredField, basis = DataBasis_pressure_w2$tBasis[,1:q3_w2])
tData_pressure_w2 <- data.frame(rbind(designW1_em[,1:4], designW2_em[,1:4]), Noise = runif(200), Coeffs_pressure_w2)

em_flow1_w2 <- BasisEmulators(tData_flow1_w2, q1_w2, mean_fn = 'step', maxdf = NULL, training_prop = 1)
par(mfrow=c(2,3), mar=c(4,4,2,2))
LeaveOneOut(em_flow1_w2[[1]]);LeaveOneOut(em_flow1_w2[[2]]);LeaveOneOut(em_flow1_w2[[3]]);LeaveOneOut(em_flow1_w2[[4]]);LeaveOneOut(em_flow1_w2[[5]])

em_flow2_w2 <- BasisEmulators(tData_flow2_w2, q2_w2, mean_fn = 'step', maxdf = NULL, training_prop = 1)
par(mfrow=c(2,3), mar=c(4,4,2,2))
LeaveOneOut(em_flow2_w2[[1]]);LeaveOneOut(em_flow2_w2[[2]]);LeaveOneOut(em_flow2_w2[[3]]);LeaveOneOut(em_flow2_w2[[4]]);LeaveOneOut(em_flow2_w2[[5]])

em_pressure_w2 <- BasisEmulators(tData_pressure_w2, q3_w2, mean_fn = 'step', maxdf = NULL, training_prop = 1)
par(mfrow=c(2,2), mar=c(4,4,2,2))
LeaveOneOut(em_pressure_w2[[1]]);LeaveOneOut(em_pressure_w2[[2]]);LeaveOneOut(em_pressure_w2[[3]])

# Predict across large LHC
Big_preds_flow1_w2 <- BasisPredGasp(BigDesign, em_flow1_w2)
Big_preds_flow2_w2 <- BasisPredGasp(BigDesign, em_flow2_w2)
Big_preds_pressure_w2 <- BasisPredGasp(BigDesign, em_pressure_w2)

# History match for each
impl_flow1_w2 <- matrix(0, nrow(BigDesign), length(l_seq))
impl_flow2_w2 <- matrix(0, nrow(BigDesign), length(l_seq))
impl_pressure_w2 <- matrix(0, nrow(BigDesign), length(l_seq))





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

# Refit basis, emulators
DataBasis_flow1_w3 <- MakeDataBasis(cbind(outputW1[,1,], outputW2[,1,], outputW3[,1,]))
q1_w3 <- ExplainT(DataBasis_flow1_w3, vtot = 0.99)
Coeffs_flow1_w3 <- CalcScores(data = DataBasis_flow1_w3$CentredField, basis = DataBasis_flow1_w3$tBasis[,1:q1_w3])
tData_flow1_w3 <- data.frame(rbind(designW1_em[,1:4], designW2_em[,1:4], designW3_em[,1:4]), Noise = runif(300), Coeffs_flow1_w3)

DataBasis_flow2_w3 <- MakeDataBasis(cbind(outputW1[,2,], outputW2[,2,], outputW3[,2,]))
q2_w3 <- ExplainT(DataBasis_flow2_w3, vtot = 0.99)
Coeffs_flow2_w3 <- CalcScores(data = DataBasis_flow2_w3$CentredField, basis = DataBasis_flow2_w3$tBasis[,1:q2_w3])
tData_flow2_w3 <- data.frame(rbind(designW1_em[,1:4], designW2_em[,1:4], designW3_em[,1:4]), Noise = runif(300), Coeffs_flow2_w3)

DataBasis_pressure_w3 <- MakeDataBasis(cbind(outputW1[,3,], outputW2[,3,], outputW3[,3,]))
q3_w3 <- ExplainT(DataBasis_pressure_w3, vtot = 0.999)
Coeffs_pressure_w3 <- CalcScores(data = DataBasis_pressure_w3$CentredField, basis = DataBasis_pressure_w3$tBasis[,1:q3_w3])
tData_pressure_w3 <- data.frame(rbind(designW1_em[,1:4], designW2_em[,1:4], designW3_em[,1:4]), Noise = runif(300), Coeffs_pressure_w3)

em_flow1_w3 <- BasisEmulators(tData_flow1_w3, q1_w3, mean_fn = 'step', maxdf = NULL, training_prop = 1)
par(mfrow=c(2,3), mar=c(4,4,2,2))
LeaveOneOut(em_flow1_w3[[1]]);LeaveOneOut(em_flow1_w3[[2]]);LeaveOneOut(em_flow1_w3[[3]]);LeaveOneOut(em_flow1_w3[[4]]);LeaveOneOut(em_flow1_w3[[5]])

em_flow2_w3 <- BasisEmulators(tData_flow2_w3, q2_w3, mean_fn = 'step', maxdf = NULL, training_prop = 1)
par(mfrow=c(2,3), mar=c(4,4,2,2))
LeaveOneOut(em_flow2_w3[[1]]);LeaveOneOut(em_flow2_w3[[2]]);LeaveOneOut(em_flow2_w3[[3]]);LeaveOneOut(em_flow2_w3[[4]]);LeaveOneOut(em_flow2_w3[[5]])

em_pressure_w3 <- BasisEmulators(tData_pressure_w3, q3_w3, mean_fn = 'step', maxdf = NULL, training_prop = 1)
par(mfrow=c(2,2), mar=c(4,4,2,2))
LeaveOneOut(em_pressure_w3[[1]]);LeaveOneOut(em_pressure_w3[[2]]);LeaveOneOut(em_pressure_w3[[3]])

# Predict across large LHC
Big_preds_flow1_w3 <- BasisPredGasp(BigDesign, em_flow1_w3)
Big_preds_flow2_w3 <- BasisPredGasp(BigDesign, em_flow2_w3)
Big_preds_pressure_w3 <- BasisPredGasp(BigDesign, em_pressure_w3)

# History match for each
impl_flow1_w3 <- matrix(0, nrow(BigDesign), 4)
impl_flow2_w3 <- matrix(0, nrow(BigDesign), 4)
impl_pressure_w3 <- matrix(0, nrow(BigDesign), 4)

for (k in 1:4){
  Big_impl_flow1_w3 <- HistoryMatch(DataBasis_flow1_w3, obs[,1] - DataBasis_flow1_w3$EnsembleMean, Big_preds_flow1_w3$Expectation, Big_preds_flow1_w3$Variance, Error = Sigma_e[,,k], Disc = 0*diag(512), weightinv = Inv[[k]])
  impl_flow1_w3[,k] <- Big_impl_flow1_w3$impl
  
  Big_impl_flow2_w3 <- HistoryMatch(DataBasis_flow2_w3, obs[,2] - DataBasis_flow2_w3$EnsembleMean, Big_preds_flow2_w3$Expectation, Big_preds_flow2_w3$Variance, Error = Sigma_e[,,k], Disc = 0*diag(512), weightinv = Inv[[k]])
  impl_flow2_w3[,k] <- Big_impl_flow2_w3$impl
  
  Big_impl_pressure_w3 <- HistoryMatch(DataBasis_pressure_w3, obs[,3] - DataBasis_pressure_w3$EnsembleMean, Big_preds_pressure_w3$Expectation, Big_preds_pressure_w3$Variance, Error = Sigma_e[,,k], Disc = 0*diag(512), weightinv = Inv[[k]])
  impl_pressure_w3[,k] <- Big_impl_pressure_w3$impl
}

apply(impl_flow1_w3 < qchisq(0.995, 512), 2, sum)
apply(impl_flow2_w3 < qchisq(0.995, 512), 2, sum)
apply(impl_pressure_w3 < qchisq(0.995, 512), 2, sum)


