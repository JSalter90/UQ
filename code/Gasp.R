source('code/basis.R')
source('code/FastHM.R')

# From ExeterUQ https://github.com/BayesExeter/ExeterUQ - edit path to read this
current_wd <- getwd()
setwd('~/Dropbox/ExeterUQ')
source('BuildEmulator/BuildEmulator.R')
setwd(current_wd)

library(RobustGaSP)

#' Building a single GaSP emulator for basis output
#' 
#' Given tData object, fit an emulator for the selected coefficient
#' 
#' @param Response a string indicating which output is being emulated. Must be consistent with a column of tData
#' @param tData data frame containing (in order): a) the design, b) a column containing noise, c) the basis coefficients
#' @param mean_fn the structure allowed in the mean function. If NULL, fits a constant mean model. If 'linear', fits a model with a linear term in each of the inputs. If 'step', selects a mean function via stepwise regression
#' @param training_prop proportion of the data to use to fit the model, sampled at random
#' @param Fouriers if fitting a mean function via step, should Fourier terms be considered?
#' @param linModel defaults to NULL. If not, gives an object of type 'lm' to be used as the mean function
#' @param nugget should a nugget be estimated? Defaults to TRUE
#' @param maxdf maximum number of terms allowed in the mean functiuon, if step used. Defaults to 0.1*size of training data
#' 
#' @return \item{em}{An rgasp emulator}
#' \item{em_lm}{If mean_fn = 'step', the regression model that was fitted}
#' \item{active}{If mean_fn = 'step', the variables that were deemed to be active}
#' \item{type}{Label indicating that the emulator was fitted with rgasp}
#' \item{mean_fn}{Label indicating the type of mean function used}
#' \item{train_data}{The subset of the data that was used to fit the emulator}
#' \item{validation_data}{The subset of the data that was not used. If training_prop = 1, this is empty}
#' 
#' @export
BuildGasp <- function(Response, tData, mean_fn = NULL, training_prop = 0.75, Fouriers = FALSE, linModel = NULL, nugget = TRUE, maxdf = NULL, ...){
  lastCand <- which(names(tData)=="Noise")
  n <- dim(tData)[1]
  ind_response <- which(names(tData) == Response)
  
  # Split into training and validation data randomly
  if (training_prop < 1){
    n_t <- ceiling(training_prop * n)
    inds_t <- sample(1:n, n_t)
    train_input <- tData[inds_t, 1:(lastCand-1)]
    train_response <- tData[inds_t, ind_response]
    validation_input <- tData[-inds_t, 1:(lastCand-1)]
    validation_response <- tData[-inds_t, ind_response]
    if (!is.null(mean_fn)){
      if (mean_fn == 'step'){
        tData_train <- tData[inds_t,]
      }
    }
  }
  
  else {
    n_t <- n
    train_input <- tData[, 1:(lastCand-1)]
    train_response <- tData[, ind_response]
    validation_input <- NULL
    validation_response <- NULL
    if (!is.null(mean_fn)){
      if (mean_fn == 'step'){
        tData_train <- tData
      }
    }
  }

  # 'Trend' option defines mean function
  # e.g. have column of 1s for intercept, thereafter parameters
  # So for linear mean, just repeat design after column of 1s
  if (is.null(mean_fn)){
    em <- rgasp(design = train_input, response = train_response, nugget.est = nugget, ...)
  }
  
  else if (mean_fn == 'linear'){
    X <- cbind(rep(1,n_t), train_input)
    em <- rgasp(design = train_input, response = train_response, trend = as.matrix(X), nugget.est = nugget, ...)
  }
  
  else if (mean_fn == 'step'){
    em_lm <- EMULATE.lm(Response, tData_train, dString="tData", tcands = names(tData)[1:lastCand], tcanfacs=NULL, maxOrder=NULL, TryFouriers = Fouriers, maxdf = maxdf)
    # Extract design matrix for this linear model
    X <- model.matrix(em_lm$linModel)
    # Active variables
    if (Fouriers == TRUE){
      active <- names(tData)[1:(lastCand-1)]
    }
    else {
      active <- em_lm$Names
    }
    if (length(active) > 0){
      train_input <- train_input[,active]
      validation_input <- validation_input[,active]
    }
    if (length(active) == 1){
      train_input <- as.matrix(train_input, ncol = 1)
      colnames(train_input) <- active
      validation_input <- as.matrix(validation_input, ncol = 1)
      colnames(validation_input) <- active
    }
    em <- rgasp(design = train_input, response = train_response, trend = as.matrix(X), nugget.est = nugget, ...)
  }
  
  else if (mean_fn == 'lm'){
    em_lm <- mean_fn
    #X <- model.matrix(linModel)
    # In case linear model wasn't fitted with exact same data:
    tt <- terms(linModel)
    tt <- delete.response(tt)
    mm <- model.frame(tt, train_input)
    X <- model.matrix(tt, mm)
    em <- rgasp(design = train_input, response = train_response, trend = as.matrix(X), nugget.est = nugget, ...)
  }
  
  train_data <- cbind(train_input, train_response)
  colnames(train_data)[dim(train_data)[2]] <- Response
  validation_data <- cbind(validation_input, validation_response)
  if (training_prop < 1){
    colnames(validation_data)[dim(train_data)[2]] <- Response
  }
  
  if (!is.null(mean_fn)){
    if (mean_fn == 'step'){
      return(list(em = em, lm = em_lm, active = active, type = 'rgasp', mean_fn = mean_fn, train_data = train_data, validation_data = validation_data))
    }
    else if (mean_fn == 'lm'){
      return(list(em = em, lm = linModel, type = 'rgasp', mean_fn = mean_fn, train_data = train_data, validation_data = validation_data))
    }
    else {
      return(list(em = em, type = 'rgasp', mean_fn = mean_fn, train_data = train_data, validation_data = validation_data))
    }
  }
  if (is.null(mean_fn)){
    return(list(em = em, type = 'rgasp', mean_fn = mean_fn, train_data = train_data, validation_data = validation_data))
  }
}

#' Evaluating GaSP emulator predictions
#' 
#' Given an object output by BuildGasp, makes predictions for a set of inputs
#' 
#' @param Design a data frame containing the input parameters, where each row is a point at which to evaluate the emulator
#' @param emulator an object output by BuildGasp (requires mean_fn, active etc.)
#' 
#' @return an object containing the mean, standard deviation, and lower and upper bounds of the 95% posterior credible interval (see predict,rgasp-method)
#' 
#' @export
PredictGasp <- function(Design, emulator){
  if (is.null(emulator$mean_fn)){
    noise_ind <- colnames(Design) == 'Noise'
    if (sum(noise_ind) > 0){
      noise_ind <- which(colnames(Design) == 'Noise')
      Design <- Design[,c(1:(noise_ind-1))]
    }
    preds <- predict(emulator$em, Design)
  }
  
  else if (emulator$mean_fn == 'linear'){
    noise_ind <- colnames(Design) == 'Noise'
    if (sum(noise_ind) > 0){
      noise_ind <- which(colnames(Design) == 'Noise')
      Design <- Design[,c(1:(noise_ind-1))]
    }
    X <- cbind(rep(1,dim(Design)[1]), Design)
    preds <- predict(emulator$em, Design, testing_trend = as.matrix(X))
  }
  
  else if (emulator$mean_fn == 'lm'){
    tt <- terms(emulator$lm)
    Terms <- delete.response(tt)
    mm <- model.frame(Terms, Design, xlev = emulator$lm$xlevels)
    X <- model.matrix(Terms, mm, contrasts.arg = emulator$lm$contrasts)
    preds <- predict(emulator$em, Design, testing_trend = as.matrix(X))
  }
  
  else if (emulator$mean_fn == 'step'){
    # Need to make sure that the new design aligns with the ordering given by active variables
    active <- emulator$active
    if (length(active) > 0){
      Design <- Design[,active]
    }
    if (length(active) == 1){
      Design <- as.matrix(Design, ncol = 1)
      colnames(Design) <- active
      Design <- as.data.frame(Design)
    }
    # Create trend matrix using linear model
    tt <- terms(emulator$lm$linModel)
    Terms <- delete.response(tt)
    mm <- model.frame(Terms, Design, xlev = emulator$lm$linModel$xlevels)
    X <- model.matrix(Terms, mm, contrasts.arg = emulator$lm$linModel$contrasts)
    # Can now predict as before
    preds <- predict(emulator$em, Design, testing_trend = as.matrix(X))
  }
  return(preds)
}


#' Fitting multiple GaSP emulators
#' 
#' A function for fitting emulators to the coefficients on each vector of the truncated basis
#' 
#' @param tData a data frame containing the design, a noise column, and projections of the field onto the basis
#' @param HowManyEmulators number of coefficients to build emulators for
#' 
#' @return A list of emulators
#' 
#' @export
BasisEmulators <- function(tData, HowManyEmulators, type = 'rgasp', ...){
  lastCand <- which(names(tData)=="Noise")
  tfirst <- lastCand + 1
  if(is.null(HowManyEmulators))
    HowManyEmulators <- length(names(tData)) - lastCand
  #if (type == 'Stan'){
  #  lapply(1:HowManyEmulators, function(k) try(BuildStanEmulator(Response=names(tData)[lastCand+k], tData=tData, cands=names(tData)[1:lastCand], 
  #                                                               additionalVariables=additionalVariables[[k]], maxdf=ceiling(length(tData[,1])/10)+1, 
  #                                                               sigmaPrior = sigmaPrior, nuggetPrior = nuggetPrior, activePrior = activePrior, 
  #                                                               activeVariables = activeVariables, prior.params = prior.params, ...), silent = TRUE))
  #}
  if (type == 'rgasp' | type == 'gasp' | type == 'Gasp'){
    lapply(1:HowManyEmulators, function(k) BuildGasp(Response = names(tData)[lastCand+k], tData = tData, ...))
  }
}

InitialBasisEmulators <- BasisEmulators

#' Predicting for several emulators simultaneously
#'
#' Given a list of emulators for basis vectors, evaluates each
#' 
#' @param Design a data frame containing the points at which to evaluate the emulators
#' @param emulators a list of BuildGasp emulators
#' 
#' @return \item{Expectation}{the posterior expectations}
#' \item{Variance}{the posterior variances}
#' 
#' @export
BasisPredGasp <- function(Design, emulators){
  require(parallel)
  EmOutput <- mclapply(1:length(emulators), function(e) PredictGasp(Design,emulators[[e]]))
  Expectation <- Variance <- matrix(0, nrow = dim(Design)[1], ncol = length(emulators))
  for (j in 1:length(emulators)){
    Expectation[,j] <- EmOutput[[j]]$mean
    Variance[,j] <- EmOutput[[j]]$sd^2
  }
  return(list(Expectation = Expectation, Variance = Variance))
}

#' Prediction and History Matching for GaSP emulators
#'
#' Takes emulators, evaluates expectations and variances for space-filling design, and history matches
#'
#' @param DataBasis object containing the basis used in emulation ($tBasis)
#' @param Obs observation vector (length \ell), must be centred
#' @param Ems a list of BuildGasp emulators
#' @param tData matrix containing parameter values
#' @param ns number of parameter settings to evaluate emulators at
#' @param Error observation error variance matrix
#' @param Disc discrepancy variance matrix
#' @param weightinv if not NULL, the inverse of W = var_err + var_disc, used for projection
#' @param Design if not NULL, passes a design at which to evaluate emulators and implausibility
#' @param PreviousWave if not NULL, provides the output of a previous PredictAndHM object, and evaluates the current NROY points from the previous design
#' #' 
#' @return \item{Design}{Space-filling design of ns points at which the emulators were evaluated}
#' \item{Expectation}{Emulator expectations}
#' \item{Variance}{Emulator variances}
#' \item{impl}{Vector of implausibilities corresponding to the rows of Expectation and Variance}
#' \item{bound}{The chi-squared bound for an \ell-dimensional field}
#' \item{nroy}{Percentage of parameter settings that are not ruled out, using bound}
#' \item{inNROY}{Vector indicating whether a parameter setting is ruled out}
#'
#' @export
PredictAndHM <- function(DataBasis, Obs, Ems, tData, ns = 1000, Error, Disc, weightinv = NULL, BasisUncertainty = FALSE, input_range = c(-1,1),
                         Design = NULL, PreviousWave = NULL){
  T_f <- qchisq(0.995, dim(DataBasis$tBasis)[1]) # bound for ruling out on field
  npar <- which(colnames(tData) == "Noise") - 1

  design_flag <- !(is.null(Design))
  
  if (is.null(Design) & is.null(PreviousWave)){
    Design <- (input_range[2] - input_range[1]) * as.data.frame(randomLHS(ns, npar)) + input_range[1]
    colnames(Design) <- colnames(tData)[1:npar]
  }
  
  if (!(is.null(PreviousWave))){
    inNROY_inds <- which(PreviousWave$inNROY == TRUE)
    Design <- PreviousWave$Design[inNROY_inds,] # only evaluate at not ruled out points
  }
  
  #if (Ems[[1]]$type == 'Stan'){
  #  EmOutput <- lapply(1:length(Ems), function(e) EMULATOR.gpstan(Design,Ems[[e]], FastVersion = TRUE))
  #  Expectation <- Variance <- matrix(0, nrow = dim(Design)[1], ncol = length(Ems))
  #  for (j in 1:length(Ems)){
  #    Expectation[,j] <- EmOutput[[j]]$Expectation
  #    Variance[,j] <- EmOutput[[j]]$Variance
  #  }
  #}
  
  if (Ems[[1]]$type == 'rgasp'){
    Preds <- BasisPredGasp(Design, Ems)
  }
  
  FieldHM <- HistoryMatch(DataBasis, Obs, Preds$Expectation, Preds$Variance, Error, Disc, weightinv = weightinv)
  
  if (!(is.null(PreviousWave))){
    FieldHM$nroy <- PreviousWave$nroy * FieldHM$nroy
  }
  
  if (design_flag == TRUE){
    print("Proportion of given design not ruled out:")
  }
  else {
    print("Proportion of original space not ruled out:")
  }
  print(FieldHM$nroy)
  return(list(Design = Design, Expectation = Preds$Expectation, Variance = Preds$Variance, impl = FieldHM$impl, bound = FieldHM$bound, nroy = FieldHM$nroy, inNROY = FieldHM$inNROY))
}

#' Prediction and history matching for multiple fields
#'
#' Given a basis, emulator, observations etc. for multiple fields, predicts and history matches for each
#'
#' @param DataBasis list of DataBasis objects
#' @param Obs list of centred observation vectors
#' @param Ems list of emulator lists
#' @param tData list of tData objects
#' @param ns size of design to sample
#' @param Error list of observation error variance matrices
#' @param Disc list of discrepancy variance matrices
#' @param weightinv if not NULL, a list of (var_err + var_disc)^{-1} for each field
#' @param Design if not NULL, passes a design at which to evaluate emulators and implausibility
#' @param PreviousWave if not NULL, provides the output of a previous PredictAndHM object, and evaluates the current NROY points from the previous design
#'
#' @return \item{Design}{Space-filling design of ns points at which the emulators were evaluated}
#' \item{Expectation}{A list of emulator expectations for each field}
#' \item{Variance}{A list of emulator variances}
#' \item{impl}{A matrix of implausibilities, with rows corresponding to Design, column corresponding to each field}
#' \item{bound}{Vector with the chi-squared bound for each field}
#' \item{nroy}{Percentage of parameter settings that are not ruled out, using bound, for each field individually}
#' \item{inNROY}{Matrix corresponding to impl and bound, indicating whether each combination of parameter setting and field is ruled out}
#'
#' @export
PredictAndHM_multi <- function(DataBasis, Obs, Ems, tData, ns = 1000, Error, Disc, weightinv = NULL, BasisUncertainty = FALSE, input_range = c(-1,1),
                               Design = NULL, PreviousWave = NULL){
  m <- length(DataBasis)
  if (!(length(Obs) == m)){
    stop('Different number of observations and DataBasis objects provided - check all lists have same length')
  }
  if (!(length(Ems) == m)){
    stop('Different number of emulators and DataBasis objects provided - check all lists have same length')
  }
  if (!(length(tData) == m)){
    stop('Different number of tData objects and DataBasis objects provided - check all lists have same length')
  }
  if (!(length(Error) == m)){
    stop('Different number of Error matrices and DataBasis objects provided - check all lists have same length')
  }
  if (!(length(Disc) == m)){
    stop('Different number of Discrepancy matrices and DataBasis objects provided - check all lists have same length')
  }
  
  if (is.null(PreviousWave)){
    output <- PredictAndHM(DataBasis[[1]], Obs[[1]], Ems[[1]], tData[[1]], ns, Error[[1]], Disc[[1]], weightinv[[1]], BasisUncertainty, input_range)
    output1 <- mclapply(2:m, function(e) PredictAndHM(DataBasis[[e]], Obs[[e]], Ems[[e]], tData[[e]], ns, Error[[e]], Disc[[e]], weightinv[[e]], BasisUncertainty, input_range, Design = output$Design, PreviousWave = PreviousWave))
    output <- c(list(output), output1)
  }
  
  if (!is.null(PreviousWave)){
    output <- mclapply(1:m, function(e) PredictAndHM(DataBasis[[e]], Obs[[e]], Ems[[e]], tData[[e]], ns, Error[[e]], Disc[[e]], weightinv[[e]], BasisUncertainty, input_range, PreviousWave = PreviousWave))
  }
  
  # Re-format the output
  Design <- output[[1]]$Design
  Expectation <- lapply(1:m, function (e) output[[e]]$Expectation)
  Variance <- lapply(1:m, function (e) output[[e]]$Variance)
  impl <- matrix(unlist(lapply(1:m, function (e) output[[e]]$impl)), nrow = dim(Design)[1], ncol = m)
  bound <- unlist(lapply(1:m, function (e) output[[e]]$bound))
  nroy <- unlist(lapply(1:m, function (e) output[[e]]$nroy))
  inNROY <- matrix(unlist(lapply(1:m, function (e) output[[e]]$inNROY)), nrow = dim(Design)[1], ncol = m)
  
  return(list(Design = Design, Expectation = Expectation, Variance = Variance, impl = impl, bound = bound, nroy = nroy, inNROY = inNROY))
}





#' Include NROY points from previous waves
#' 
#' At wave k > 1, add any points from previous waves that are in the current NROY space, and define the new DataBasis object and centred observations by adding these to the new ensemble
#' 
#' @param DesignHM output from running PredictAndHM over the old design
#' @param NewDesign new inputs
#' @param NewData output corresponding to NewDesign
#' @param DataBasis object from previous wave
#' @param Obs centred observation from previous wave
#' 
#' @return \item{Design}{Design with new wave, followed by any NROY runs from previous waves}
#' \item{DataBasis}{Basis, centred data, and ensemble mean for the new ensemble combined with previous NROY runs}
#' \item{Obs}{Centred observations, given new ensemble mean} 
#' 
#' @export
AddPreviousNROY <- function(DesignHM, NewDesign, NewData, DataBasis, Obs, ...){
  NROYinds <- which(DesignHM$inNROY == TRUE)
  FullDesign <- rbind(NewDesign, DesignHM$Design[NROYinds,])
  OldRawData <- DataBasis$CentredField[,NROYinds] + DataBasis$EnsembleMean
  RawObs <- Obs + DataBasis$EnsembleMean
  AllData <- cbind(NewData, OldRawData)
  NewDataBasis <- MakeDataBasis(data = AllData, RemoveMean = TRUE, ...)
  NewObs <- RawObs - NewDataBasis$EnsembleMean
  return(list(Design = FullDesign, DataBasis = NewDataBasis, Obs = NewObs))
}




