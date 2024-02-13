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
