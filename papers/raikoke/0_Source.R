#### Loading general required for NAME/Raikoke 2019 emulation ####
library(ggplot2)
library(dplyr)
library(reshape2)
library(fields)
library(viridis)
library(plyr)
library(cowplot)
library(GGally)
library(akima)

source('code/mogp_save.R')
source('code/mogp_new_wrappers.R')
source('code/mogp_basis.R')

# From ExeterMOGP
source('BuildEmulator/BuildEmulator.R')

# From ExeterUQ



# Defining parameter ranges
minMER <- 50.7 * 10^7 * (9 - 0.551)^(1/0.241) * 0.33 * 0.5/100
maxMER <- 50.7 * 10^7 * (17 - 0.551)^(1/0.241) * 3 * 20/100
parRanges <- data.frame(parameter = colnames(design)[3:13],
                        lower = c(9,0.5,0.33,minMER,1350,0,9,0.0025,100,0.27,log(minMER)),
                        upper = c(17,20,3,maxMER,2500,17,15,2.75,900,1.74,log(maxMER)))



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


# History matching experiment, treating each left-out run as observations in turn
PseudoExperiment <- function(tData, val_inds, em_pred, obs_error){
  
  n <- length(val_inds)
  val_data <- tData[val_inds,]
  
  overall_impl <- numeric(n) # implausibility of assumed obs in each experiment, for overall emulator
  total_matches <- numeric(n) # how many METs are considered not implausible, in each experiment
  overall_size <- pseudo_size <- cons_size <- numeric(n) # size of the different spaces
  
  for (i in 1:n){
    pseudo_obs <- val_data$LogTotal[i]
    
    # Calculate overall implausibility
    overall_impl[i] <- abs(em_pred$overall$Mean[val_inds[i]] - pseudo_obs) / sqrt(em_pred$overall$SD[val_inds[i]]^2 + obs_error)
    
    # Size of overall NROY space
    impl_all <- (abs(em_pred$overall$Mean - pseudo_obs) / sqrt(em_pred$overall$SD^2 + obs_error))[val_inds]
    overall_size[i] <- sum(impl_all < 3) / n 
    
    # Calculate implausibility for each MET
    impl_MET <- matrix(0, n, 18)
    for (j in 1:18){
      impl_MET[,j] <- abs(em_pred$met[[j]]$Mean[val_inds] - pseudo_obs) / sqrt(em_pred$met[[j]]$SD[val_inds]^2 + obs_error)
    }
    
    total_matches[i] <- sum(impl_MET[i,] < bound) # count how many METs the chosen 'obs' are in NROY for 
    
    # Size of pseudo NROY space
    pseudo_size[i] <- sum(apply(impl_MET < 3, 1, sum) > 8) / n
    
    # Size of conservative NROY space
    cons_size[i] <- sum(apply(impl_MET < 3, 1, sum) > 0) / n
  }
  
  return(list(overall_impl,
              overall_size,
              total_matches,
              pseudo_size,
              cons_size))
}



