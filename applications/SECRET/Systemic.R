# SECRET competition, September 2023
# Code for systemic model, Exeter team

# Exact results may differ when run this code, as some randomness involved - 
# selecting initial training dataset
# fitting mean function uses a noise vector (shouldn't make much/any difference)
# sampling wave 2 from current NROY
# sampling new points at wave 3/4

# We've reproduced the exact same results up to end of wave 2

# See Example_EmulateTS.html for more detail on the emulation approach, applied to a single time series (generalises to other high dimensional output, e.g., stacking together multiple series as here)

library(R.matlab)
library(reshape2)

# Basis emulation functions are found at https://github.com/JSalter90/UQ
# Also reads a file from https://github.com/BayesExeter/ExeterUQ, edit paths in Gasp.R to reflect location of this
setwd('~/Dropbox/UQ/') # edit directory on local machine to where https://github.com/JSalter90/UQ is cloned
source('code/Gasp.R')

# Load inputs
setwd('~/Dropbox/UQ/applications/SECRET')
design_sim <- readRDS('data/systemic/design_sim.rds') # on original scale
design_em <- readRDS('data/systemic/design_em.rds') # on [-1,1]^5

# Load data
flow_all <- readRDS("data/systemic/flow_all.rds") # 7 flow series across design, 512x1250x7
maxmin_all <- readRDS("data/systemic/maxmin_all.rds") # max/min pressure across design, 1250x3

# Stacking flows to emulate simultaneously (ignored the max/min pressure, but these could be added to the end of the output vector)
flow_stacked <- aperm(flow_all, c(1,3,2)) # changing to 512x7x1250
dim(flow_stacked) <- c(512*7, nrow(design_em)) # vectorising output

# Load observations
obs <- c(readMat('data/systemic/obs.mat')[[1]]) # also vectorising, ignoring pressure

# At the time, the true inputs were unknown
# Included here so can plot truth alongside history matching results
truth <- data.frame(f2 = -32.9,
                      f3 = 426000,
                      fs2 = -40.6,
                      fs3 = 643000,
                      alpha = 0.88)

# Scaled
truth_em <- truth
truth_em[,1] <- (truth_em[,1] + 45) / ((20)/2) - 1
truth_em[,2] <- (truth_em[,2] - 2*10^5) / ((7*10^5)/2) - 1
truth_em[,3] <- (truth_em[,3] + 45) / ((20)/2) - 1
truth_em[,4] <- (truth_em[,4] - 2*10^5) / ((7*10^5)/2) - 1
truth_em[,5] <- (truth_em[,5] - 0.85) / ((0.09)/2) - 1

# We know the 'true' implausibility of all runs (i.e, we use the model output, with emulator variance = 0)
true_impl <- numeric(ncol(flow_stacked))
for (i in 1:length(true_impl)){
  true_impl[i] <- sum((obs - flow_stacked[,i])^2) # assuming obs error = 1
}
inNROY <- which(true_impl < qchisq(0.995, 512*7))
design_em[inNROY,]

# Split into train/validation
n <- 200 # number of training points
set.seed(37281)
samp <- sample(1:nrow(design_em), nrow(design_em))
train_inds <- samp[1:n]
val_inds <- samp[-c(1:n)]

train_design <- design_em[train_inds,]
val_design <- design_em[val_inds,]

train_full <- flow_stacked[,train_inds]
val_full <- flow_stacked[,val_inds]

# Construct basis, project, create emulator data
DataBasis_full <- MakeDataBasis(train_full) # centres the data, calculates basis
q_full <- ExplainT(DataBasis_full, vtot = 0.99) # vectors required to explain 99% (arbitrary; chosen here after experimentation showed that it was possible to emulate these vectors, and that adding more added little)
Coeffs_full <- Project(data = DataBasis_full$CentredField, basis = DataBasis_full$tBasis[,1:q_full]) # project centred data onto basis
tData_full <- data.frame(train_design[,1:5], Noise = runif(n), Coeffs_full) # combining inputs with coefficients

# Fit emulator
# As above, these settings chosen based on prior experimentation with the systemic model
em_full <- BasisEmulators(tData_full, q_full, mean_fn = 'step', maxdf = NULL, training_prop = 1)

# LOO cross-validation, could also consider predicting on validation set
par(mfrow=c(2,2), mar=c(4,4,2,2))
LeaveOneOut(em_full[[1]]);LeaveOneOut(em_full[[2]]);LeaveOneOut(em_full[[3]]);LeaveOneOut(em_full[[4]])
LeaveOneOut(em_full[[5]]);LeaveOneOut(em_full[[6]]);LeaveOneOut(em_full[[7]]);LeaveOneOut(em_full[[8]])

# Broadly ok, nothing systematically wrong, and generally ok close to the observations:
obs_coeffs <- Project(obs - DataBasis_full$EnsembleMean,
                      DataBasis_full$tBasis[,1:q_full])

# Predict across large LHC
# BigDesign <- 2*as.data.frame(randomLHS(100000, 5)) - 1
# colnames(BigDesign) <- colnames(design_em)[1:5]
# Load from file for consistency
BigDesign <- readRDS('data/systemic/BigDesign.rds')

# Predict across design
Big_preds1 <- BasisPredGasp(BigDesign, em_full)

# History match
# $Expectation and $Variance should be matrices with size (number of rows in BigDesign) x q
# j^th column should be Expectation (Variance) of j^th coefficient
# This code efficiently calculates implausibility (over original field), but to do so requires argument weightinv = W^{-1}, where W = Var_e + Var_{disc}
# To allow this to be efficient, weightinv must have a certain structure (it requires attributes flagging whether it is identity, diagonal or more complex)
# However, if you generate this inverse using GetInverse, it will automatically generate this in the correct form
Err <- 1*diag(512*7) # uncorrelated error, with arbitrary variance as in theory this is zero here
Winv <- GetInverse(Err) # creating inverse
Big_impl1 <- HistoryMatch(DataBasis_full, 
                          obs - DataBasis_full$EnsembleMean, # compare relative to ensemble mean, as this is removed prior to emulation
                          Big_preds1$Expectation, 
                          Big_preds1$Variance, 
                          Error = Err, 
                          Disc = 0*Err, # discrepancy is zero, provides zero matrix with correct dimension
                          weightinv = Winv)
Big_impl1$nroy # proportion of input space in NROY
Big_impl1$bound # taken from chi-squared with ell = 512*7 degrees of freedom
summary(Big_impl1$impl)

# Visualising NROY
library(GGally)
BigDesign$NROY <- Big_impl1$inNROY
k <- 1:5
p <- ggpairs(BigDesign[1:5000,], columns=k,
             ggplot2::aes(color=NROY) , upper = list(continuous = wrap("density", alpha = 0.5), combo = "box_no_facet"),
             lower = list(continuous = wrap("points", alpha = 0.3), combo = wrap("dot_no_facet", alpha = 0.4)),
             diag = list(continuous = wrap("densityDiag", alpha = 0.3)),
             legend = 1) +
  theme(legend.position = "bottom") + scale_colour_manual(values = c("#F8766D", "#00BFC4")) +
  scale_fill_manual(values = c("#F8766D", "#00BFC4"))
# Add truth
for(i in 1:length(k)) {
  p1 <- getPlot(p, i, i) + geom_vline(xintercept = as.numeric(truth_em[k[i]]))
  p <- putPlot(p, p1, i, i)
}
p

# Instead, evaluate at known runs
Val_preds1 <- BasisPredGasp(val_design, em_full)
Val_impl1 <- HistoryMatch(DataBasis_full, obs - DataBasis_full$EnsembleMean, Val_preds1$Expectation, Val_preds1$Variance, Error = Err, Disc = 0*diag(dim(Err)[1]), weightinv = Winv)
Val_impl1$nroy

# Minimum implausibility across 100k samples, known runs, vs truth
BigDesign[which.min(Big_impl1$impl),]
val_design[which.min(Val_impl1$impl),]
truth_em

# Unsurprisingly, not perfect after only 1 wave - emulator variance not negligible - only used 200 training points
# However, now we can sample from NROY, and train new emulators with these new runs
# The emulators should be more accurate, as now have denser training samples in the region of parameter space closer to the observations
# We've also likely removed parts of space that lead to very different behaviour, hence no longer trying to capture this with the emulators

#### Wave 2 ####
# Ordinarily we'd want to sample from NROY generally, and run new sim,ulations
# Because of the time-limited nature of SECRET, and Matlab license/internet connection issues, we instead sampled the 200 new points from the existing validation set
inNROY1 <- which(Val_impl1$inNROY)

nroy_design <- val_design[inNROY1,]
nroy_flows <- val_full[,inNROY1]

set.seed(581) # a seed was not used during the 3 hour timeframe of the competition
w2_inds <- sample(1:nrow(nroy_design), n)

train_full_w2 <- nroy_flows[,w2_inds]
val_full_w2 <- nroy_flows[,-w2_inds]

train_design_w2 <- nroy_design[w2_inds,]
val_design_w2 <- nroy_design[-w2_inds,]

# New basis
DataBasis_full_w2 <- MakeDataBasis(train_full_w2)
q_full_w2 <- ExplainT(DataBasis_full_w2, vtot = 0.99)
Coeffs_full_w2 <- Project(data = DataBasis_full_w2$CentredField, basis = DataBasis_full_w2$tBasis[,1:q_full_w2])
tData_full_w2 <- data.frame(train_design_w2[,1:5], Noise = runif(n), Coeffs_full_w2)

# As a sanity check - does our basis better represent the observations now?
ReconError(obs - DataBasis_full$EnsembleMean, basis = DataBasis_full$tBasis[,1:q_full], scale = FALSE)
ReconError(obs - DataBasis_full_w2$EnsembleMean, basis = DataBasis_full_w2$tBasis[,1:q_full_w2], scale = FALSE)

# New emulators
em_full_w2 <- BasisEmulators(tData_full_w2, q_full_w2, mean_fn = 'step', maxdf = NULL, training_prop = 1)
par(mfrow=c(3,3), mar=c(4,4,2,2))
LeaveOneOut(em_full_w2[[1]]);LeaveOneOut(em_full_w2[[2]]);LeaveOneOut(em_full_w2[[3]]);LeaveOneOut(em_full_w2[[4]])
LeaveOneOut(em_full_w2[[5]]);LeaveOneOut(em_full_w2[[6]]);LeaveOneOut(em_full_w2[[7]]);LeaveOneOut(em_full_w2[[8]]);LeaveOneOut(em_full_w2[[9]])

# Predict across large LHC
Big_preds1_w2 <- BasisPredGasp(BigDesign, em_full_w2)

# History match
# The choice of 'observation error' is again arbitrary here as really it's zero
Err_w2 <- 0.1*diag(512*7)
Winv_w2 <- GetInverse(Err_w2)
Big_impl1_w2 <- HistoryMatch(DataBasis_full_w2, 
                             obs - DataBasis_full_w2$EnsembleMean, 
                             Big_preds1_w2$Expectation, 
                             Big_preds1_w2$Variance, 
                             Error = Err_w2, 
                             Disc = 0*Err_w2, 
                             weightinv = Winv_w2)
# Need points in NROY at both w1 and w2:
#sum(Big_impl1$inNROY & Big_impl1_w2$inNROY) # actual amount in this space is arbitrary, because using arbitrary variance. Only care about smallest values for this application

# Plot
# BigDesign$NROY2 <- Big_impl1$inNROY & Big_impl1_w2$inNROY
# p <- ggpairs(BigDesign[1:5000,], columns=k,
#              ggplot2::aes(color=NROY2) , upper = list(continuous = wrap("density", alpha = 0.5), combo = "box_no_facet"),
#              lower = list(continuous = wrap("points", alpha = 0.3), combo = wrap("dot_no_facet", alpha = 0.4)),
#              diag = list(continuous = wrap("densityDiag", alpha = 0.3)),
#              legend = 1) +
#   theme(legend.position = "bottom") + scale_colour_manual(values = c("#F8766D", "#00BFC4")) +
#   scale_fill_manual(values = c("#F8766D", "#00BFC4"))
# # Add truth
# for(i in 1:length(k)) {
#   p1 <- getPlot(p, i, i) + geom_vline(xintercept = as.numeric(truth_em[k[i]]))
#   p <- putPlot(p, p1, i, i)
# }
# p

# Points in input space that minimise impl at w2, vs truth
BigDesign[which.min(Big_impl1_w2$impl),] # and in fact input 76468 is the closest it is possible to get to the truth from the 100k samples in BigDesign
order(Big_impl1_w2$impl)[1:10]
Big_impl1$impl[order(Big_impl1_w2$impl)[1:10]] # these 'best' 10 all in W1 NROY as well
BigDesign[order(Big_impl1_w2$impl)[1:10],]
truth_em

# Would usually now sample from NROY, do new simulations, refit emulators, etc. (as emulator uncertainty can likely still be reduced, so should be able to zoom in further)
# Due to issues with running new simulations (University Matlab licensing issues) and time constraints, only performed 10 new simulations
# These were chosen as the 10 points in the 100k of BigDesign that minimised the W2 implausibility
# Ordinarily, would do space-filling in NROY, but aim of competition is to find run as close as possible to z

#### Wave 3 ####
# The 10 selected points from above may not be exactly the same, however should generally identify same region of input space
# Here, load in the 10 new simulations
flow_W3 <- readRDS('data/systemic/flow_W3.rds')
design_W3 <- readRDS('data/systemic/design_W3.rds')

# Have we identified better simulations at each wave?
true_impl_w1 <- numeric(n)
for (i in 1:length(true_impl_w1)){
  true_impl_w1[i] <- sum((obs - train_full[,i])^2)
}

true_impl_w2 <- numeric(n)
for (i in 1:length(true_impl_w2)){
  true_impl_w2[i] <- sum((obs - train_full_w2[,i])^2)
}

true_impl_w3 <- numeric(nrow(design_W3))
for (i in 1:length(true_impl_w3)){
  true_impl_w3[i] <- sum((obs - flow_W3[,i])^2)
}

# Most summaries decrease - on average finding simulations closer to the truth
summary(true_impl_w1)
summary(true_impl_w2) # the 'best' (i.e. minimum impl) simulation in wave 2 is not better than wave 1. Ordinarily we wouldn't expect this, however we were sampling from a finite set of 1250 runs, and had by chance already included the run with minimum impl out of these 1250. Hence W2 design was guaranteed to not reduce this  
summary(true_impl_w3)


#### Wave 4 ####
# Designed by some small perturbations around best point from wave 3, due to lack of time
flow_W4 <- readRDS('data/systemic/flow_W4.rds')
design_W4 <- readRDS('data/systemic/design_W4.rds')

true_impl_w4 <- numeric(nrow(design_W4))
for (i in 1:length(true_impl_w4)){
  true_impl_w4[i] <- sum((obs - flow_W4[,i])^2)
}

# Generally closer again, as well as finding run with lower implausibility
summary(true_impl_w4)
