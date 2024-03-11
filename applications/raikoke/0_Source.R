#### Loading general required for NAME/Raikoke 2019 emulation ####
packages <- c('ggplot2','dplyr','cowplot','reticulate','pracma', 'invgamma', 'GenSA', 'far', 'fields', 'lhs', 'maps', 
              'mco', 'mvtnorm', 'ncdf4', 'parallel', 'reshape2', 'shape', 'tensor', 'viridis', 'withr', 'loo','MASS')
sapply(packages, require, character.only = TRUE, quietly = TRUE)

source('code/mogp_save.R')
source('code/mogp_new_wrappers.R')
source('code/mogp_basis.R')

# Emulation here done with mogp_emulator, found here: https://github.com/alan-turing-institute/mogp_emulator
mogp_dir <- '/Users/jamessalter/mogp_emulator/' # location of mogp installation
# This section is handled within BuildEmulator.R (below), given mogp_dir
# twd <- getwd()
# setwd(mogp_dir)
# mogp_emulator <- import("mogp_emulator")
# mogp_priors <- import("mogp_emulator.Priors")
# mogp_kernels <- import("mogp_emulator.Kernel")

# We use parts of the the R front end from https://bayesexeter.github.io/ExeterUQ_MOGP/
exeterUQ_mogp_dir <- '/Users/jamessalter/Documents/ExeterUQ_MOGP'
uqdir <- getwd()
setwd(exeterUQ_mogp_dir)
source('BuildEmulator/BuildEmulator.R')
setwd(uqdir)

# Load design matrix
design <- readRDS('applications/raikoke/data/design.rds')

# Defining parameter ranges
minMER <- 50.7 * 10^7 * (9 - 0.551)^(1/0.241) * 0.33 * 0.5/100
maxMER <- 50.7 * 10^7 * (17 - 0.551)^(1/0.241) * 3 * 20/100
parRanges <- data.frame(parameter = colnames(design)[3:13],
                        lower = c(9,0.5,0.33,minMER,1350,0,9,0.0025,100,0.27,log(minMER)),
                        upper = c(17,20,3,maxMER,2500,17,15,2.75,900,1.74,log(maxMER)))

# Scaling design to [-1,1]
ScaledDesign <- design[,c(3:7,9:13)] # ignore MET
for (i in 1:ncol(ScaledDesign)){
  ind <- which(parRanges$parameter == colnames(ScaledDesign)[i])
  ScaledDesign[,i] <- ScaledDesign[,i] - parRanges$lower[ind]
  ScaledDesign[,i] <- ScaledDesign[,i] / ((parRanges$upper[ind] - parRanges$lower[ind])/2)
  ScaledDesign[,i] <- ScaledDesign[,i] - 1
}
ScaledDesign$MER <- NULL # for consistency; never use this variable in emulation/calibration as log(MER) is included

# Observation information
# For history matching, require estimated obs (mean/median), observation error variance
# $Mean gives mean (on log scale), $Var gives estimated variance (on log scale)
# Just load this in
obs <- readRDS('applications/raikoke/data/obs.rds')

#### Additional functions for this work specifically ####
# Plot NAME output
PlotPlume <- function(output, xlims = c(140,200), ylims = c(40, 61), zlims = c(0,log(5000)), longitude = lon, latitude = lat, obs = NULL, cols = viridis(100), ...){
  lon_ind <- which(between(longitude, xlims[1], xlims[2]))
  lat_ind <- which(between(latitude, ylims[1], ylims[2]))
  full_grid <- expand.grid(lon = lon_ind, lat = lat_ind)
  full_grid <- full_grid %>% left_join(output[,c('xind', 'yind', 'output')], by = c('lon' = 'xind', 'lat' = 'yind'))
  image.plot(sort(lon[lon_ind]), sort(lat[lat_ind]), matrix(log(full_grid$output+1), length(lon_ind), length(lat_ind))[order(lon[lon_ind]),], 
             col = cols, zlim = zlims, xlim = xlims, ylim = ylims, ...)
  map('world2', add = TRUE)
  points(source_loc[1], source_loc[2], col = 'red', pch = 17, cex = 2)
  if (!is.null(obs)){
    points(0.5*(obs$X_Lon1+obs$X_Lon2), 0.5*(obs$Y_Lat1+obs$Y_Lat2), pch = 4, col = 'yellow', cex = 0.25)
  }
}

# Validate both the overall and MET-specific emulators
ValidateBoth <- function(em1, em2, tData, Design, TrainInds, ValInds, TrainIndsMET, ValIndsMET){
  
  TrainData <- tData[TrainInds,]
  ValData <- tData[ValInds,]
  
  val_preds_total <- ValidationMOGP(NewData = ValData, em1, tData = TrainData, which.emulator = 1, 
                                    ParamNames = colnames(TrainData)[em1$fitting.elements$ActiveIndices[[1]]])
  val_preds_total <- as.data.frame(val_preds_total)
  colnames(val_preds_total) <- c('Mean', 'Lower', 'Upper')
  val_preds_total$Truth <- ValData[,ncol(ValData)]
  val_preds_total$In95 <- val_preds_total$Truth >= val_preds_total$Lower & val_preds_total$Truth <= val_preds_total$Upper
  val_preds_total$MET <- Design$MET[ValInds]
  print(sum(val_preds_total$In95) / length(val_preds_total$In95))
  
  val_preds_total_MET <- NULL
  for (m in 0:17){
    inds <- which(Design$MET[ValIndsMET] == m)
    tmp <- ValidationMOGP(NewData = tData[ValIndsMET[which(Design$MET[ValIndsMET] == m)],], em2[[m+1]], tData = tData[TrainIndsMET[which(Design$MET[TrainIndsMET] == m)],], which.emulator = 1, 
                          ParamNames = colnames(tData[TrainIndsMET,])[em2[[m+1]]$fitting.elements$ActiveIndices[[1]]])
    tmp <- as.data.frame(tmp)
    colnames(tmp) <- c('Mean', 'Lower', 'Upper')
    tmp$Truth <- tData[ValIndsMET[which(Design$MET[ValIndsMET] == m)], ncol(tData)]
    tmp$MET <- m
    val_preds_total_MET <- rbind(val_preds_total_MET, tmp)
  }
  
  val_preds_total_MET$In95 <- val_preds_total_MET$Truth >= val_preds_total_MET$Lower & val_preds_total_MET$Truth <= val_preds_total_MET$Upper
  print(sum(val_preds_total_MET$In95) / length(val_preds_total_MET$In95))
  
  tmp <- rbind(data.frame(val_preds_total, Type = 'Overall'),
               data.frame(val_preds_total_MET, Type = 'MET-specific'))

  plot1 <- ggplot(tmp, aes(x = Truth, y = Mean, col = as.factor(MET))) + 
    geom_errorbar(aes(ymin = Lower, ymax = Upper)) +
    geom_point() +
    facet_wrap(vars(Type)) +
    geom_abline(slope = 1) +
    labs(col = 'm', y = 'Emulator', x = 'NAME output')
  
  return(list(plot1 = plot1,
              overall = val_preds_total,
              met = val_preds_total_MET))
}

# Predict for both the overall and MET-specific emulators
PredictBoth <- function(em1, em2, Design){
  
  Pred1 <- em1$mogp$predict(as.matrix(Design))
  Pred1 <- data.frame(Design,
                      Mean = Pred1$mean[1,],
                      SD = sqrt(Pred1$unc[1,]))
  Pred1$Lower <- Pred1$Mean - 1.96*Pred1$SD
  Pred1$Upper <- Pred1$Mean + 1.96*Pred1$SD
  
  Pred_tmp <- NULL
  for (i in 1:18){
    Pred_tmp[[i]] <- em2[[i]]$mogp$predict(as.matrix(Design))
  }
  
  Pred2 <- NULL
  for (i in 1:18){
    Pred2[[i]] <- data.frame(Design,
                             Mean = Pred_tmp[[i]]$mean[1,],
                             SD = sqrt(Pred_tmp[[i]]$unc[1,]),
                             Lower = Pred_tmp[[i]]$mean[1,] - 1.96*sqrt(Pred_tmp[[i]]$unc[1,]),
                             Upper = Pred_tmp[[i]]$mean[1,] + 1.96*sqrt(Pred_tmp[[i]]$unc[1,]),
                             MET = i-1)
  }
  
  return(list(overall = Pred1,
              met = Pred2))
}


#' History matching experiment, treating each left-out run as observations in turn
#'
#' Works for single or multiple outputs, loops over all validation points
#'
#' @param tData design matrix, including output in column `LogTotal`. Can be a list for multiple outputs.
#' @param val_inds vector corresponding to left-out points when fitting emulator. Can be a list.
#' @param em_pred emulator predictions (containing `$Mean`, `$SD`) corresponding to `tData`. Can be a list.
#' @param obs_error vector of observation error variances corresponding to each output
#' @param kmax kth max implausibility, defaults to 1 (i.e. take max implausibility across outputs)
#'
#' @return \item{overall_impl}{Overall implausibility for each of the left-out runs}
#' \item{overall_size}{Size of overall NROY (across the left-out runs) in each experiment}
#' \item{total_matches}{For each left-out run in turn, number of METs for which this point is in NROY}
#' \item{pseudo_size}{Size of pseudo NROY in each experiment}
#' \item{cons_size}{Size of conservative NROY in each experiment}
#'
#' @export
PseudoExperiment <- function(tData, val_inds, em_pred, obs_error, kmax = 1, bound = 3){
  
  if (length(dim(tData)) == 0){
    ell <- length(tData)
    stopifnot(ell == length(em_pred)) # checking provided same number of sets of predictions as number of outputs
    stopifnot(ell == length(val_inds)) # checking provided validation indices for each output
    stopifnot(ell == length(obs_error)) # checking provided obs error for each output
  }
  
  else {
    ell <- 1
  }
  
  if (ell == 1){
    tData <- list(tData) # converting to list so can use common code below for all ell
    val_inds <- list(val_inds)
    em_pred <- list(em_pred)
  }
  
  n <- length(val_inds[[1]])
  
  val_data <- list()
  for (k in 1:ell){
    val_data[[k]] <- tData[[k]][val_inds[[k]],]
  }

  overall_impl <- numeric(n) # implausibility of assumed obs in each experiment, for overall emulator
  total_matches <- numeric(n) # how many METs are considered not implausible, in each experiment
  overall_size <- pseudo_size <- cons_size <- numeric(n) # size of the different spaces
  
  for (i in 1:n){
    pseudo_obs <- numeric(ell)
    for (k in 1:ell){
      pseudo_obs[k] <- val_data[[k]]$LogTotal[i]
    }

    # Calculate overall implausibility
    tmp_impl <- numeric(ell)
    for (k in 1:ell){
      tmp_impl[k] <- abs(em_pred[[k]]$overall$Mean[val_inds[[k]][i]] - pseudo_obs[k]) / sqrt(em_pred[[k]]$overall$SD[val_inds[[k]][i]]^2 + obs_error[k])
    }
    
    overall_impl[i] <- kth_max(tmp_impl, k = kmax)
    
    # Size of overall NROY space
    impl_all <- matrix(NA, n, ell)
    for (k in 1:ell){
      impl_all[,k] <- (abs(em_pred[[k]]$overall$Mean - pseudo_obs[k]) / sqrt(em_pred[[k]]$overall$SD^2 + obs_error[k]))[val_inds[[k]]]
    }
    impl_all <- apply(impl_all, 1, kth_max, k = kmax)
    overall_size[i] <- sum(impl_all < bound) / n
    
    # Calculate implausibility for each MET
    impl_MET <- array(0, dim = c(n, 18, ell))
    for (k in 1:ell){
      for (j in 1:18){
        impl_MET[,j,k] <- abs(em_pred[[k]]$met[[j]]$Mean[val_inds[[k]]] - pseudo_obs[k]) / sqrt(em_pred[[k]]$met[[j]]$SD[val_inds[[k]]]^2 + obs_error[k])
      }
    }
    impl_MET <- apply(impl_MET, c(1,2), kth_max, k = kmax)
    total_matches[i] <- sum(impl_MET[i,] < bound) # count how many METs the chosen 'obs' are in NROY for 
    
    # Size of pseudo NROY space
    pseudo_size[i] <- sum(apply(impl_MET < bound, 1, sum) > 8) / n
    
    # Size of conservative NROY space
    cons_size[i] <- sum(apply(impl_MET < bound, 1, sum) > 0) / n
  }
  
  return(list(overall_impl = overall_impl,
              overall_size = overall_size,
              total_matches = total_matches,
              pseudo_size = pseudo_size,
              cons_size = cons_size))
}

SummariseExperiment <- function(PseudoExperiment){
  data.frame(Type = c('Pseudo', 'Overall', 'Cons'),
             Errors = c(sum(PseudoExperiment$total_matches < 9),
                        sum(PseudoExperiment$overall_impl > 3), 
                        sum(PseudoExperiment$total_matches == 0)),
             Size = c(median(PseudoExperiment$pseudo_size), 
                      median(PseudoExperiment$overall_size),
                      median(PseudoExperiment$cons_size)))
}

kth_max <- function(x,k) {
  sorted_values <- sort(x, decreasing = TRUE)
  sorted_values[k]  # returns the kth maximum
}


