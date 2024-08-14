library(ggplot2)
library(PLNmodels)
library(dplyr)
library(viridis)
library(hetGP)
library(reshape2)

#' Finding the Poisson Log Normal PCA basis for a matrix of counts, corresponding to n simulations with ell outputs
#'
#' @param data an ell x n matrix of counts, with each row corresponding to an output location, and each column corresponding to model output
#' @param ... other arguments to pass to `PLNPCA`
#'
#' @return \item{tBasis}{The latent basis, an orthogonal basis of dimension ell x rank, given by the SVD over the (centred) latent model outputs}
#' \item{Data}{A copy of the data used in the CountBasis call}
#' \item{EnsembleMean}{The mu parameter in the PLNPCA decomposition, which is not necessarily the mean of the data}
#' \item{LatentMean}{The mean of BB^T in the PLNPCA decomposition (there is no orthonormality constraint in the estimation step, hence prior to finding tBasis, need to remove this latent mean)}
#' \item{Coeffs}{The coefficients for the n runs on the latent basis}
#' \item{PLN}{Storing the PLNPCA object, containing all estimated parameters}
#'
#' @examples
#' @export
CountBasis <- function(data, ...){
  t1 <- Sys.time()
  n <- ncol(data)
  ell <- nrow(data)
  data_fit <- PLNmodels::prepare_data(as.data.frame(t(data)), rep(NA, n), 'none')
  pln_basis <- PLNmodels::PLNPCA(Abundance ~ 1, data  = data_fit, ...)
  t2 <- Sys.time()
  print(t2 - t1)

  best <- PLNmodels::getBestModel(pln_basis)

  # Format in standard way, but also store raw outputs
  EnsembleMean <- as.numeric(coef(best)) # not quite the ensemble mean, but similar
  Coeffs <- best$scores
  tBasis <- best$rotation
  LatentMean <- as.numeric(apply(best$latent_pos, 2, mean)) # the latent part (relative to ens mean) is not necessarily mean zero
  Data <- data

  # Relabel coefficient columns
  colnames(Coeffs)[1:ncol(Coeffs)] <- paste("C",1:ncol(Coeffs),sep="")

  return(list(tBasis = tBasis,
              Data = Data,
              EnsembleMean = EnsembleMean,
              LatentMean = LatentMean,
              Coeffs = Coeffs,
              PLN = pln_basis))
}



#' Joining inputs and outputs
#'
#' Formats the input and outputs for later use in emulation
#'
#' @param Design n x p matrix or data frame, where each row corresponds to the p input parameters for an individual simulation
#' @param CountBasis output provided by CountBasis. Uses the $Coeffs output, the ordering here must be consistent with the provided Design
#' @param q How many vectors to include coefficients from. Defaults to ncol(CountBasis$Coeffs)
#' @param Noise whether to include a noise vector (used in emulation). Defaults to TRUE
#'
#' @return Inputs and outputs required for emulation
#'
#' @export
GetEmDataCount <- function(Design, CountBasis, q = NULL, Noise = TRUE){
  if (Noise){
    Noise <- runif(nrow(Design), -1, 1)
  }
  if (is.null(q)){
    q <- ncol(CountBasis$Coeffs)
  }
  tData <- cbind(Design, Noise, CountBasis$Coeffs[,1:q])
  return(tData)
}


#' Fitting multiple HetGP emulators
#'
#' A function for fitting emulators to the coefficients on each vector of the truncated basis
#'
#' @param tData a data frame containing the design, a noise column, and projections of the field onto the basis
#' @param HowManyEmulators number of coefficients to build emulators for
#' @param ... other arguments to pass to `BuildHet`
#'
#' @return A list of emulators
#'
#' @export
BasisEmulatorsHet <- function(tData, HowManyEmulators, ...){
  lastCand <- which(names(tData)=="Noise")
  tfirst <- lastCand + 1
  if(is.null(HowManyEmulators)){
    HowManyEmulators <- length(names(tData)) - lastCand
  }
  lapply(1:HowManyEmulators, function(k) BuildHet(Response = names(tData)[lastCand+k], tData = tData, ...))
}


#' Building a single hetGP emulator for basis output
#'
#' Given tData object, fit an emulator for the selected coefficient
#'
#' @param Response a string indicating which output is being emulated. Must be consistent with a column of tData
#' @param tData data frame containing (in order): a) the design, b) a column containing noise, c) the basis coefficients
#' @param training_prop proportion of the data to use to fit the model, sampled at random
#' @param ... other arguments to pass to `mleHetGP`
#'
#' @return \item{em}{An HetGP emulator}
#' \item{type}{Label indicating that the emulator was fitted with HetGP}
#' \item{train_data}{The subset of the data that was used to fit the emulator}
#' \item{validation_data}{The subset of the data that was not used. If training_prop = 1, this is empty}
#'
#' @export
BuildHet <- function(Response, tData, training_prop = 0.75, ...){
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
  }

  else {
    n_t <- n
    train_input <- tData[, 1:(lastCand-1)]
    train_response <- tData[,ind_response]
    validation_input <- NULL
    validation_response <- NULL
  }

  # For hetGP, need to process tData to handle replicates
  het_input <- hetGP::find_reps(X = as.matrix(train_input), Z = train_response)

  # Fit emulator
  em <- hetGP::mleHetGP(X = list(X0 = het_input$X0, Z0 = het_input$Z0, mult = het_input$mult),
                        Z = het_input$Z, covtype = 'Matern3_2', ...)

  train_data <- cbind(train_input, train_response)
  colnames(train_data)[dim(train_data)[2]] <- Response
  validation_data <- cbind(validation_input, validation_response)
  if (training_prop < 1){
    colnames(validation_data)[dim(train_data)[2]] <- Response
  }

  return(list(em = em,
              type = 'het',
              train_data = train_data,
              validation_data = validation_data))
}




#' Evaluating hetGP emulator predictions
#'
#' Given an object output by BuildHet, makes predictions for a set of inputs
#'
#' @param Design a data frame containing the input parameters, where each row is a point at which to evaluate the emulator
#' @param emulator an object output by BuildHet
#'
#' @return an object containing the mean, variance (with nugget removed), and nugget variance, at each input location
#'
#' @export
PredictHet <- function(Design, emulator){
  
  # Ensure columns are ordered in the same way as when the emulator was trained
  col_names <- colnames(emulator$train_data)[-ncol(emulator$train_data)]
  Design <- Design[,col_names]
  
  preds <- predict(emulator$em, as.matrix(Design))
  return(preds)
}


#' Validating hetGP
#'
#' @param emulator an object output by BuildHet
#' @param ValidationData locations to predict at
#' @param IndivPars Create plots for each input
#'
#' @return an object containing the mean, variance (with nugget removed), and nugget variance, at each input location
#'
#' @export
ValidateHet <- function(emulator, ValidationData = NULL, IndivPars = FALSE){

  if (!is.null(ValidationData)){
    emulator$validation_data <- ValidationData[,colnames(emulator$train_data)]
  }
  resp_ind <- dim(emulator$validation_data)[2]
  design <- emulator$validation_data[,-resp_ind]
  if (length(emulator$active) == 1){
    design <- as.matrix(design, ncol = 1)
    colnames(design) <- emulator$active
  }
  response <- emulator$validation_data[,resp_ind]
  preds <- PredictHet(design, emulator)

  vars <- preds$sd2 + preds$nugs
  preds$lower95 <- preds$mean - 1.96*sqrt(vars)
  preds$upper95 <- preds$mean + 1.96*sqrt(vars)
  upp <- max(c(preds$upper95, response))
  low <- min(c(preds$lower95, response))
  preds$truth <- response

  preds$sd2var <- preds$cov <- NULL

  preds$In95 <- preds$truth >= preds$lower95 & preds$truth <= preds$upper95
  perc_outside <- round(sum(preds$In95 == FALSE) / length(preds$In95) * 100, 1)
  cols <- c('darkgrey', viridis::viridis(100)[31], viridis::viridis(100)[81])

  # Ensuring good points still coloured green if no points outside
  if (perc_outside == 0){
    cols[2:3] <- viridis::viridis(100)[81]
  }

  plots <- ggplot(as.data.frame(preds), aes(x = .data$truth, y = .data$mean, col = .data$In95)) +
    geom_errorbar(aes(ymin = .data$lower95, ymax = .data$upper95), col = cols[1]) +
    geom_point() +
    scale_colour_manual(values = c(cols[2:3])) +
    geom_abline(slope = 1, alpha = 0.6) +
    labs(y = 'Prediction', x = 'Truth', title = paste0('Outside 95% = ', perc_outside, '%')) +
    theme(legend.position = 'none')

  if (IndivPars == TRUE){
    plots_inputs <- NULL
    for (j in 1:dim(design)[2]){
      plot_data <- data.frame(x = design[,j], as.data.frame(preds))
      plots_inputs[[j]] <- ggplot(plot_data, aes(x = .data$x, y = .data$mean, col = .data$In95)) +
        geom_errorbar(aes(ymin = .data$lower95, ymax = .data$upper95), col = cols[1]) +
        geom_point() +
        scale_colour_manual(values = c(cols[2:3])) +
        labs(y = 'Prediction', x = paste0(colnames(design)[j])) +
        theme(legend.position = 'none')
    }
  }

  if (IndivPars == FALSE){
    return(plots)
  }

  else{
    return(list(plot = plots,
                input = plots_inputs))
  }
}





#' Predicting for several emulators simultaneously
#'
#' Given a list of emulators for basis vectors, evaluates each
#'
#' @param Design a data frame containing the points at which to evaluate the emulators
#' @param emulators a list of BuildHet emulators
#'
#' @return \item{Expectation}{the posterior expectations}
#' \item{Variance}{the posterior variances}
#'
#' @export
BasisPredHet <- function(Design, emulators){
  EmOutput <- parallel::mclapply(1:length(emulators), function(e) PredictHet(Design,emulators[[e]]))
  Expectation <- Variance <- matrix(0, nrow = dim(Design)[1], ncol = length(emulators))
  for (j in 1:length(emulators)){
    Expectation[,j] <- EmOutput[[j]]$mean
    Variance[,j] <- EmOutput[[j]]$sd2 + EmOutput[[j]]$nugs
  }
  return(list(Expectation = Expectation, Variance = Variance))
}


#' Basis reconstruction of training data
#'
#' For original matrix of counts, reconstructs the field given a basis. Does the equivalent of fitted() or predict() for PLNmodels, but by hand
#'
#' @param CountBasis output provided by CountBasis
#' @param q The number of basis vectors to use, defaults to NULL (uses all vectors)
#' @param inds Which simulations to reconstruct, by giving a vector of row indices. Defaults to NULL (uses all rows)
#'
#' @return Reconstruction for each simulation
#'
#' @export
ReconFitted <- function(CountBasis, q = NULL, inds = NULL){
  ell <- nrow(CountBasis$tBasis)

  if (is.null(q)){
    q <- ncol(CountBasis$tBasis)
  }

  if (is.null(inds)){
    n <- nrow(CountBasis$Coeffs)
    inds <- 1:n
  }

  else {
    n <- length(inds)
  }

  C <- CountBasis$PLN$models[[1]]$model_par$C[,1:q]
  S2 <- CountBasis$PLN$models[[1]]$var_par$S2[inds,1:q]

  var <- matrix(0, ell, n)

  if (q > 1){
    for (i in 1:n){
      var[,i] <- diag(C %*% diag(S2[i,]) %*% t(C))
    }
  }

  else {
    for (i in 1:n){
      var[,i] <- diag(S2[i] * C %*% t(C))
    }
  }

  latent <- CountBasis$EnsembleMean + CountBasis$LatentMean + CountBasis$tBasis[,1:q] %*% t(CountBasis$Coeffs[inds,1:q])
  response <- exp(latent + 0.5*var)

  return(response)
}


#' Basis reconstruction from coefficients
#'
#' Given a vector of coefficients (and their variances) on the latent basis, reconstructs the original field
#' Mostly used within other functions and applied over a matrix of expectations and variances, i.e. from a prediction
#'
#' @param Coeffs a vector of coefficients, usually the expectation at an input x
#' @param Variance the associated variance on the coefficients
#' @param CountBasis output of `CountBasis()`, assuming that the coefficients relate to the 1st q vectors of the $tBasis field
#'
#' @return Reconstructed expectation over the original field
#'
#' @export
ReconExpectation <- function(Coeffs, Variance, CountBasis){
  q <- length(Coeffs)
  Basis <- CountBasis$tBasis[,1:q]

  if (q == 1){
    var <- diag(Variance * Basis %*% t(Basis))
    latent <- CountBasis$EnsembleMean + CountBasis$LatentMean + Basis * c(Coeffs)
  }
  else {
    var <- diag(Basis %*% diag(Variance) %*% t(Basis))
    latent <- CountBasis$EnsembleMean + CountBasis$LatentMean + Basis %*% c(Coeffs)
  }

  response <- exp(latent + 0.5*var)
  return(c(response))
}

#' Basis reconstruction from predictions
#'
#' Given expectation/variance across set of inputs, calculates expectation of each on original field
#'
#' @param BasisPred Expectation (in $Expectation) and variance (in $Variance) at prediction locations. Both must be matrices of dimension \eqn{n_t \times q}, where \eqn{n_t} is the number of prediction locations, and \eqn{q} is the number of emulated basis vectors.
#' @param CountBasis Output of `CountBasis()`, containing the basis used for emulation/projection/reconstruction
#' @param q Number of basis vectors to use. By default, will use all found in BasisPred.
#'
#' @return Reconstructed expectations over the original field, for each input
#'
#' @export
ReconPred <- function(BasisPred, CountBasis, q = NULL){
  ell <- nrow(CountBasis$tBasis)
  ns <- nrow(BasisPred$Expectation)
  # If haven't set q, use all emulated vectors to reconstruct
  if (is.null(q)){
    q <- ncol(BasisPred$Expectation)
  }
  # If q is still NULL, there was only 1 vector in Predictions
  if (is.null(q)){
    q <- 1
  }
  mu <- BasisPred$Expectation
  var <- BasisPred$Variance
  if (q == 1){
    recons <- lapply(1:ns, function(i) ReconExpectation(mu[i], var[i], CountBasis))
  }
  else {
    recons <- lapply(1:ns, function(i) ReconExpectation(mu[i,1:q], var[i,1:q], CountBasis))
  }
  recons <- matrix(unlist(recons), ell)
  return(recons)
}



#' Sampling from a PLNPCA basis emulator
#'
#' Given the mean and variance on the latent basis at a set of locations, draws samples of the reconstructed, original field
#'
#' @param BasisPred Expectation (in $Expectation) and variance (in $Variance) at prediction locations. Both must be matrices of dimension \eqn{n_t \times q}, where \eqn{n_t} is the number of prediction locations, and \eqn{q} is the number of emulated basis vectors.
#' @param CountBasis Output of `CountBasis()`, containing the basis used for emulation/projection/reconstruction
#' @param ns Number of samples to draw at each input, defaults to 1000
#' @param ReturnAll (logical) Should all individual samples be returned, or only the mean and 95% interval? Defaults to FALSE to save memory.
#' @param BasisUncertainty Include uncertainty from deleted basis vectors when sampling, defaults to TRUE
#'
#' @return \item{mean}{A matrix of dimension \eqn{\ell \times n_t}, where each row contains the mean across the `ns` samples of the original field for a given input set}
#' \item{lower}{The 2.5th percentile across the samples at each input set}
#' \item{upper}{The 97.5th percentile across the samples at each input set}
#' \item{samples}{If `ReturnAll = TRUE`, an \eqn{\ell \times n_s \times n_t} array containing the individual samples}
#'
#' @export
CountBasisEmSamples <- function(BasisPred, CountBasis, ns = 1000, ReturnAll = FALSE, BasisUncertainty = TRUE){
  n <- nrow(BasisPred$Expectation)
  q <- ncol(BasisPred$Expectation)
  if (is.null(n)){ # i.e. a single vector has been provided
    n <- 1
    q <- length(BasisPred$Expectation)
  }

  ell <- dim(CountBasis$tBasis)[1]
  Basis <- CountBasis$tBasis[,1:q]

  if (BasisUncertainty){
    EstVar <- apply(CountBasis$Coeffs, 2, var)[-c(1:q)] # variance on these vectors
    
    # Want full basis when reconstructing now
    Basis <- CountBasis$tBasis
    
    # Append zero means and these variances to $Expectation, $Variance
    if (n == 1){
      BasisPred$Expectation <- c(BasisPred$Expectation, rep(0, ncol(Basis)-q))
      BasisPred$Variance <- c(BasisPred$Variance, EstVar)
    }
    
    if (n > 1){
      BasisPred$Expectation <- cbind(BasisPred$Expectation, matrix(0, n, ncol(Basis) - q))
      BasisPred$Variance <- cbind(BasisPred$Variance, matrix(rep(EstVar, each = n), n, ncol(Basis) - q))
    }
    
    q <- ncol(Basis)
  }

  if (n == 1){
    samp <- matrix(rnorm(q*ns,
                         mean = rep(BasisPred$Expectation),
                         sd = rep(sqrt(BasisPred$Variance))), nrow = ns, byrow = TRUE)
    rec <- CountBasis$EnsembleMean + CountBasis$LatentMean + Basis %*% t(samp)
    em_samp <- matrix(rpois(ell*ns, c(exp(rec))), nrow = ell)
  }

  if (n > 1){
    em_samp <- array(0, dim = c(ell, ns, n))
    for (i in 1:n){
      samp <- matrix(rnorm(q*ns,
                           mean = rep(BasisPred$Expectation[i,]),
                           sd = rep(sqrt(BasisPred$Variance[i,]))), nrow = ns, byrow = TRUE)
      rec <- CountBasis$EnsembleMean + CountBasis$LatentMean + Basis %*% t(samp)
      em_samp[,,i] <- matrix(rpois(ell*ns, c(exp(rec))), nrow = ell)
    }
  }

  if (n == 1){
    samp_mean <- apply(em_samp, 1, mean)
    lower <- apply(em_samp, 1, quantile, probs = 0.025)
    upper <- apply(em_samp, 1, quantile, probs = 0.975)
  }

  if (n > 1){
    samp_mean <- apply(em_samp, c(1,3), mean)
    lower <- apply(em_samp, c(1,3), quantile, probs = 0.025)
    upper <- apply(em_samp, c(1,3), quantile, probs = 0.975)
  }

  if (ReturnAll){
    return(list(mean = samp_mean,
                lower = lower,
                upper = upper,
                samples = em_samp))
  }

  else {
    return(list(mean = samp_mean,
                lower = lower,
                upper = upper))
  }
}




#' Aggregating samples of the original field
#'
#' Often want to aggregate each sample to a total, or across a set of locations. Takes samples drawn from the emulator and reconstructed to the original field, and aggregates and calculates summaries of these samples.
#'
#' @param Samples A set of samples from the PLNPCA emulator, reconstructed to the original field. Usually the `$samples` field of `CountBasisEmSamples()`
#' @param locs Indices of locations (across the output field) to aggregate over. Defaults to NULL, which uses all outputs
#'
#' @return \item{mean}{Mean across the `ns` samples of the original field, aggregated to the given locations, for each input set}
#' \item{lower}{The 2.5th percentile across the (aggregated) samples at each input set}
#' \item{upper}{The 97.5th percentile across the (aggregated) samples at each input set}
#' \item{samples}{Aggregations for each sample}
#'
#' @export
AggregateSamples <- function(Samples, locs = NULL){

  ell <- dim(Samples)[1]
  ns <- dim(Samples)[2]
  # Input will either be ell x ns (i.e. for a single input), or ell x ns x n
  if (is.na(dim(Samples)[3])){
    n <- 1
  }

  else {
    n <- dim(Samples)[3]
  }

  if (is.null(locs)){
    locs <- 1:ell
  }

  # Subset to required locations, and sum
  if (length(locs) == 1){
    if (n == 1){
      sub_total <- Samples[locs,]
      lower <- quantile(sub_total, probs = 0.025)
      mu <- mean(sub_total)
      upper <- quantile(sub_total, probs = 0.975)
    }
    else {
      sub_total <- Samples[locs,,]
      lower <- apply(sub_total, 2, quantile, probs = 0.025)
      mu <- apply(sub_total, 2, mean)
      upper <- apply(sub_total, 2, quantile, probs = 0.975)
    }
  }

  else {
    if (n == 1){
      sub_data <- Samples[locs,]
      sub_total <- apply(sub_data, 2, sum)
      lower <- quantile(sub_total, probs = 0.025)
      mu <- mean(sub_total)
      upper <- quantile(sub_total, probs = 0.975)
    }

    else {
      sub_data <- Samples[locs,,]
      sub_total <- apply(sub_data, c(2,3), sum)
      lower <- apply(sub_total, 2, quantile, probs = 0.025)
      mu <- apply(sub_total, 2, mean)
      upper <- apply(sub_total, 2, quantile, probs = 0.975)
    }
  }

  return(list(mean = mu,
              lower = lower,
              upper = upper,
              samples = sub_total))
}


#' Validating predictions based on aggregations
#'
#' For a set of samples of reconstructions, and chosen output locations, compares the emulator expectation and variance across this aggregation to the true values
#'
#' @param Samples A set of samples from the PLNPCA emulator, reconstructed to the original field. Usually the `$samples` field of `CountBasisEmSamples()`
#' @param Truth The corresponding true model outputs, for comparison. Must have an ordering of rows and columns consistent with `Samples` (in terms of input and output locations)
#' @param locs Indices of locations (across the output field) to aggregate over. Defaults to NULL, which uses all outputs
#'
#' @return Validation plot comparing truth and emulator samples
#'
#' @export
ValidateSum <- function(Samples, Truth, locs = NULL){
  ell <- dim(Samples)[1]
  ns <- dim(Samples)[2]
  if (is.null(locs)){
    locs <- 1:ell
  }

  EmSum <- AggregateSamples(Samples, locs)

  if (!(nrow(Truth) == length(locs))){
    Truth <- Truth[locs,]
  }
  
  if (length(locs) == 1){
    TruthSum <- Truth
  }
  
  else {
    TruthSum <- apply(Truth, 2, sum)
  }
  
  plot_data <- data.frame(Mean = EmSum$mean,
                          Lower = EmSum$lower,
                          Upper = EmSum$upper,
                          Truth = TruthSum)

  plot_data$In95 <- plot_data$Truth >= plot_data$Lower & plot_data$Truth <= plot_data$Upper
  perc_outside <- round(sum(plot_data$In95 == FALSE) / length(plot_data$In95) * 100, 1)
  cols <- c('darkgrey', viridis::viridis(100)[31], viridis::viridis(100)[81])

  # Ensuring good points still coloured green if no points outside
  if (perc_outside == 0){
    cols[2:3] <- viridis::viridis(100)[81]
  }
  
  # Plotting zeros more clearly
  plot_data$Mean[plot_data$Mean <= 0.1] <- 0.1
  plot_data$Lower[plot_data$Lower == 0] <- 0.1
  plot_data$Upper[plot_data$Upper == 0] <- 0.1
  plot_data$Truth[plot_data$Truth == 0] <- 0.1
  
  plot <- ggplot(plot_data, aes(.data$Truth, .data$Mean, col = .data$In95)) +
    geom_errorbar(aes(ymin = .data$Lower, ymax = .data$Upper), col = cols[1]) +
    geom_point() +
    scale_colour_manual(values = c(cols[2:3])) +
    geom_abline(slope = 1, alpha = 0.6) +
    scale_x_log10(breaks = c(0.1,1,10,100,1000,10000,100000,1000000), labels = c('0','1','10','100','1000','10000','100000','1000000')) +
    scale_y_log10(breaks = c(0.1,1,10,100,1000,10000,100000,1000000), labels = c('0','1','10','100','1000','10000','100000','1000000')) +
    labs(y = 'Prediction', title = paste0('Outside 95% = ', perc_outside, '%')) +
    theme_bw() +
    theme(legend.position = 'none')

  return(plot)
}



#' 
#' 
#' Training single emulator, for running in parallel for multiple training sets
#' 
#' @param design
#' @param data 
#' @param n 
#' @param seed 
#' @param basis 
#' @param hetGP 
#' @param reps 
#' @param ... 
#' 
#' @return 
#' 
#' @export
TrainEmulator <- function(path, design, data, n, seed, id = "LAD19CD", basis = FALSE, reg = 'All', hetGP = TRUE, standardGP = FALSE, reps = FALSE, r = 10, q = 3, ...){
  
  # Create folder to write predictions/samples to
  dir.create(path)
  
  train_prop <- n / nrow(design)
  set.seed(seed)
  
  # Split replicated and unreplicated inputs in same proportion
  design_single <- subset(design, repeats == 1)
  design_reps <- subset(design, repeats > 1)
  train_inds_single <- sample(1:nrow(design_single), train_prop*nrow(design_single))
  train_inds_reps <- sample(1:nrow(design_reps), train_prop*nrow(design_reps))
  train_labels <- c(design_single$output[train_inds_single], design_reps$output[train_inds_reps])
  
  train_data <- subset(data, output %in% train_labels)
  val_data <- subset(data, !(output %in% train_labels))
  
  if (!(reps)){
    train_data <- subset(train_data, replicate == 1)
    val_data <- subset(val_data, replicate == 1)
  }
  
  # Find basis
  if (basis == 'PLNPCA'){
    # Convert data into l x n matrix
    train_data <- ProcessData(train_data, by_id = id)
    val_data <- ProcessData(val_data, by_id = id)
    
    # Find corresponding inputs
    train_design <- design[match(train_data$run$output, design$output),]
    val_design <- design[match(val_data$run$output, design$output),]
    
    train_basis <- CountBasis(train_data$data, rank = r)
    tData <- GetEmDataCount(train_design[,1:15], train_basis)
  }
  
  else if (basis == 'SVD'){
    # Convert data into l x n matrix
    train_data <- ProcessData(train_data, by_id = id)
    val_data <- ProcessData(val_data, by_id = id)
    
    # Find corresponding inputs
    train_design <- design[match(train_data$run$output, design$output),]
    val_design <- design[match(val_data$run$output, design$output),]
    
    train_basis <- MakeDataBasis(log(train_data$data+1))
    tData <- GetEmData(Design = train_design[,1:15], q = q, DataBasis = train_basis)
  }
  
  # If basis not specified, instead aggregate data and emulate overall/regional total
  else {
    train_basis <- NULL
    # Aggregate to total/region
    if (reg == 'All'){
      train_data_total <- aggregate(cbind(cumH, cumHD, cumCD, deaths) ~ output + replicate, train_data, sum)
      val_data_total <- aggregate(cbind(cumH, cumHD, cumCD, deaths) ~ output + replicate, val_data, sum)
      train_data_total$LogDeaths <- log(train_data_total$deaths)
      val_data_total$LogDeaths <- log(val_data_total$deaths)
    }
    else {
      # If given specific region
      train_data_total <- aggregate(cbind(cumH, cumHD, cumCD, deaths) ~ output + replicate, subset(train_data, region == reg), sum)
      val_data_total <- aggregate(cbind(cumH, cumHD, cumCD, deaths) ~ output + replicate, subset(val_data, region == reg), sum)
      
      # Remove any -Inf
      if (any(train_data_total$deaths == 0) | any(val_data_total$deaths == 0)){
        train_data_total$LogDeaths <- log(train_data_total$deaths + 1)
        val_data_total$LogDeaths <- log(val_data_total$deaths + 1)
      }
      
      else {
        train_data_total$LogDeaths <- log(train_data_total$deaths)
        val_data_total$LogDeaths <- log(val_data_total$deaths)
      }
    }
    
    train_design <- design[match(train_data_total$output, design$output),1:15]
    val_design <- design[match(val_data_total$output, design$output),1:15]
    tData_total <- data.frame(train_design,
                              Noise = runif(nrow(train_design), -1, 1), 
                              LogDeaths = train_data_total$LogDeaths)
    
    # Set number of emulators to be fitted equal to 1
    q <- 1
  }
  
  if (q > 1 & basis == 'PLNPCA'){
    if (standardGP){
      EmCoeffs <- BasisEmulators(tData, q, mean_fn = 'step', training_prop = 1)
      Preds_val <- BasisPredGasp(val_design[,1:15], EmCoeffs)
      val_response <- CountBasisEmSamples(Preds_val, train_basis, ns = 100, ReturnAll = TRUE)
      val_response$Run <- val_design$output
      val_response$mean <- val_response$lower <- val_response$upper <- NULL
      saveRDS(val_response, file = paste0(path, '/basis_', basis, '_em_',  seed, '.rds'))
    }
    
    if (hetGP){
      EmCoeffs <- BasisEmulatorsHet(tData, q, training_prop = 1)
      Preds_val <- BasisPredHet(val_design[,1:15], EmCoeffs)
      val_response <- CountBasisEmSamples(Preds_val, train_basis, ns = 100, ReturnAll = TRUE)
      val_response$Run <- val_design$output
      val_response$mean <- val_response$lower <- val_response$upper <- NULL
      saveRDS(val_response, file = paste0(path, '/basis_', basis, '_em_het_',  seed, '.rds'))
    }
  }
  
  if (q > 1 & basis == 'SVD'){
    if (standardGP){
      EmCoeffs <- BasisEmulators(tData, q, mean_fn = 'step', training_prop = 1)
      Preds_val <- BasisPredGasp(val_design[,1:15], EmCoeffs)
      tmp <- BasisEmSamples(Preds_val, train_basis, ns = 100, ReturnAll = TRUE)
      val_response <- NULL
      val_response$samples <- exp(tmp)-1
      val_response$samples[val_response$samples < 0] <- 0
      val_response$Run <- val_design$output
      saveRDS(val_response, file = paste0(path, '/basis_', basis, '_em_',  seed, '.rds'))
    }
    
    if (hetGP){
      EmCoeffs <- BasisEmulatorsHet(tData, q, training_prop = 1)
      Preds_val <- BasisPredHet(val_design[,1:15], EmCoeffs)
      tmp <- BasisEmSamples(Preds_val, train_basis, ns = 100, ReturnAll = TRUE)
      val_response <- NULL
      val_response$samples <- exp(tmp)-1
      val_response$samples[val_response$samples < 0] <- 0
      val_response$Run <- val_design$output
      saveRDS(val_response, file = paste0(path, '/basis_', basis, '_em_het_',  seed, '.rds'))
    }
  }
  
  if (q == 1){
    if (standardGP){
      EmTotal <- BuildGasp('LogDeaths', mean_fn = 'step', tData_total, training_prop = 1)
      Preds_val_total <- PredictGasp(val_design[,1:15], EmTotal)
      output_data <- data.frame(Run = val_data_total$output,
                                Truth = val_data_total$LogDeaths,
                                Mean = Preds_val_total$mean,
                                Lower = Preds_val_total$lower95,
                                Upper = Preds_val_total$upper95)
      output_data$In95 <- output_data$Truth >= output_data$Lower & output_data$Truth <= output_data$Upper
      saveRDS(output_data, file = paste0(path, '/', reg, '_em_', seed, '.rds'))
    }
    
    if (hetGP){
      EmTotal <- BuildHet('LogDeaths', tData_total, training_prop = 1)
      Preds_val_total <- PredictHet(val_design[,1:15], EmTotal)
      tmp_var <- Preds_val_total$sd2 + Preds_val_total$nugs
      Preds_val_total$lower <- Preds_val_total$mean-qnorm(0.975)*sqrt(tmp_var)
      Preds_val_total$upper <- Preds_val_total$mean+qnorm(0.975)*sqrt(tmp_var)
      output_data <- data.frame(Run = val_data_total$output,
                                Truth = val_data_total$LogDeaths,
                                Mean = Preds_val_total$mean,
                                Lower = Preds_val_total$lower,
                                Upper = Preds_val_total$upper)
      output_data$In95 <- output_data$Truth >= output_data$Lower & output_data$Truth <= output_data$Upper
      saveRDS(output_data, file = paste0(path, '/', reg, '_em_het_', seed, '.rds'))
    }
  }
}



# design needs to be design_em


#' For given set of emulator assumptions, trains emulators across different training/validation splits
#'
#' @param experiments 
#' @param design design, already scaled to emulator level (e.g., between -1,1)
#' @param data 
#' @param path 
#' @param ... 
#'
#' @return
#' @export
FitMulti <- function(path, experiments, design, data, id = 'LAD19CD', basis = NULL, reg = 'All', hetGP = TRUE, standardGP = FALSE, reps = FALSE, r = 10, q = 3){
  # Go through rows of experiments, contains n, seed
  parallel::mclapply(1:nrow(experiments), function(k) TrainEmulator(path, design, data, n = experiments$n[k], experiments$seed[k], id, basis, reg, hetGP, standardGP, reps, r, q),
                     mc.cores = parallel::detectCores())
}







#' Plot the latent mean
#'
#' Plots the estimated mean on the latent level (not quite equivalent to ensemble mean). Count output will generally be spatially indexed, so by default this function takes lon/lat or similar 2d coordinates
#'
#' @param CountBasis Output of `CountBasis()`, containing the basis used for emulation/projection/reconstruction
#' @param coords If provided, lon/lat locations of the outputs
#'
#' @return A plot of the latent mean
#' @import ggplot2
#' @export
PlotMean <- function(CountBasis, coords = NULL){
  # Function assumes that provided coordinates are consistent with ordering in basis
  if (!(is.null(coords))){
    plot_data <- data.frame(Longitude = coords[,1],
                            Latitude = coords[,2],
                            Mean = CountBasis$EnsembleMean)
  }

  plot <- ggplot(plot_data, aes(.data$Longitude, .data$Latitude, col = .data$Mean)) +
    geom_point(size = 3) +
    viridis::scale_colour_viridis()

  return(plot)
}

#' Plot the latent basis
#'
#' Plots the 1st q latent basis vectors. Count output will generally be spatially indexed, so by default this function takes lon/lat or similar 2d coordinates
#'
#' @param CountBasis Output of `CountBasis()`, containing the basis used for emulation/projection/reconstruction
#' @param q Number of vectors to plot
#' @param coords If provided, lon/lat locations of the outputs
#'
#' @return A plot with each panel showing a basis vector
#'
#' @export
PlotBasis <- function(CountBasis, q = 9, coords = NULL){
  ell <- nrow(CountBasis$tBasis)
  if (!(is.null(coords))){
    plot_data <- data.frame(Longitude = coords[,1],
                            Latitude = coords[,2],
                            Weight = c(CountBasis$tBasis[,1:q]),
                            Vector = rep(1:q, each = ell))
  }
  
  point_size <- ifelse(ell >= 1000, 0.25, 0.75)

  plot <- ggplot(plot_data, aes(.data$Longitude, .data$Latitude, col = .data$Weight)) +
    geom_point(size = point_size) +
    facet_wrap(vars(.data$Vector)) +
    viridis::scale_colour_viridis() +
    theme_bw()

  return(plot)
}


#' Plots reconstructions of runs from basis representation
#'
#' Plots reconstructions of runs from basis representation
#'
#' @param CountBasis Output of `CountBasis()`, containing the basis used for emulation/projection/reconstruction
#' @param q Number of latent basis vectors to use
#' @param inds Input rows to use
#' @param coords If provided, lon/lat locations of the outputs
#'
#' @return A plot of the difference between truth and basis reconstruction, for each run
#'
#' @export
PlotReconCount <- function(CountBasis, q = 1, inds = 1:16, coords = NULL){

  ell <- nrow(CountBasis$tBasis)
  basis <- CountBasis$tBasis[,1:q]
  recons <- ReconFitted(CountBasis, q, inds)
  truth <- CountBasis$Data[,inds]

  if (!(is.null(coords))){
    plot_data <- data.frame(Longitude = coords[,1],
                            Latitude = coords[,2],
                            Residual = c(truth - recons),
                            Run = rep(inds, each = ell))
  }

  lim <- max(abs(plot_data$Residual))

  plot <- ggplot(plot_data, aes(.data$Longitude, .data$Latitude, col = .data$Residual)) +
    geom_point() +
    facet_wrap(vars(.data$Run)) +
    scale_colour_gradient2(low = 'blue', mid = 'white', high = 'red', limits = c(-lim, lim))

  return(plot)
}

#' Plot coefficients on the latent basis
#'
#' Visualise variability across replicates on a given vector
#'
#' @param CountBasis Output of `CountBasis()`, containing the basis used for emulation/projection/reconstruction
#' @param k Which basis coefficients to plot
#' @param run_ids data frame with columns (output, replicate), with same ordering as CountBasis$data, CountBasis$Coeffs
#' @param plot_inds a vector for which runs to plot (i.e. if lots of runs with replicates, may want to subset these for visualiation purposes)
#'
#' @return A box plot showing variability across replicates
#'
#' @export
PlotCoeffs <- function(CountBasis, k = 1, run_ids, plot_inds = NULL){

  # In case not provided with these column names
  colnames(run_ids) <- c('output', 'replicate')

  # Extract ids of replicated runs
  rep_ids <- unique(run_ids$output[which(run_ids$replicate > 1)])

  if (is.null(plot_inds)){
    plot_inds <- 1:length(rep_ids)
  }
  
  # Extract which rows correspond to replicates of these inputs
  inds <- which(run_ids$output %in% rep_ids[plot_inds])

  # Combine required coefficients with 'output'
  plot_data <- data.frame(Coefficient = CountBasis$Coeffs[inds,k],
                          Run = run_ids$output[inds])

  plot <- ggplot(plot_data, aes(x = .data$Coefficient, y = as.factor(.data$Run))) +
    geom_boxplot(fill = 'grey') +
    labs(x = paste0('c', k), y = '') +
    geom_vline(xintercept = 0, linetype = 'dotted') +
    theme_bw() +
    theme(legend.position = "none")

  return(plot)
}


#' Plot PlNPCA emulator samples
#'
#' Plot emulator samples, vs truth if available
#'
#' @param Samples A set of samples from the PLNPCA emulator, reconstructed to the original field. Usually the `$samples` field of `CountBasisEmSamples()`
#' @param locs Indices of locations (across the output field) to aggregate over. Defaults to NULL, which uses all outputs
#' @param runs Which runs to plot
#' @param Truth The true output, if available
#'
#' @return A plot of reconstructed emulator samples, for the given runs/locations, with the true output overlaid if provided
#'
#' @export
PlotSamplesCount <- function(Samples, locs = NULL, runs = NULL, Truth = NULL){
  ell <- dim(Samples)[1]
  ns <- dim(Samples)[2]
  # Input will either be ell x ns (i.e. for a single input), or ell x ns x n
  if (is.na(dim(Samples)[3])){
    n <- 1
  }
  else {
    n <- dim(Samples)[3]
  }

  # If don't provide a subset to plot, use all
  if (is.null(locs)){
    locs <- 1:ell
  }

  # If don't provide which runs to plot, use all
  if (is.null(runs)){
    runs <- 1:n
  }

  # Subset to required locations
  if (n == 1){
    Samples <- Samples[locs,]
  }

  else {
    Samples <- Samples[locs,,runs]
  }

  plot_data <- data.frame(Input = 1:length(locs),
                          Output = c(Samples),
                          s = rep(1:ns, each = length(locs)),
                          Run = rep(runs, each = length(locs)*ns))

  line_size <- ifelse(length(runs) >= 9 | length(locs) > 50, 0.75, 1)
  
  plot_data$Output[plot_data$Output == 0] <- 0.1
  
  plot <- ggplot(plot_data, aes(.data$Input, .data$Output, col = as.factor(.data$s))) +
    geom_line() +
    geom_line(data = data.frame(aggregate(Output ~ Input + Run, plot_data, mean)), col = 'red', size = line_size) +
    geom_line(data = data.frame(aggregate(Output ~ Input + Run, plot_data, quantile, probs = 0.025)), col = 'red', size = line_size, linetype = 'dashed') +
    geom_line(data = data.frame(aggregate(Output ~ Input + Run, plot_data, quantile, probs = 0.975)), col = 'red', size = line_size, linetype = 'dashed') +
    scale_colour_manual(values = rep('grey', max(plot_data$s))) +
    facet_wrap(vars(.data$Run)) +
    scale_y_log10(breaks = c(0.1,1,10,100,1000,10000,100000), labels = c('0','1','10','100','1000','10000','100000')) +
    theme_bw() +
    theme(legend.position = 'none')
  
  if (length(locs)<=15){
    plot <- plot + scale_x_continuous(breaks=1:length(locs))
  }

  # If provided with true output to overlay
  if (!(is.null(Truth))){
    truth_data <- data.frame(Input = 1:length(locs),
                             Output = c(Truth[locs,runs]),
                             Run = rep(runs, each = length(locs)))
    truth_data$Output[truth_data$Output == 0] <- 0.1
    plot <- plot +
      geom_line(data = truth_data, col = 'black', size = line_size)
  }

  return(plot)
}


#' Plot replicated training input
#'
#' For a single input, plots sampled fields and all replicates
#'
#' @param Samples A set of samples from the PLNPCA emulator, reconstructed to the original field. Usually the `$samples` field of `CountBasisEmSamples()`
#' @param Replicates The true replicates
#' @param locs Indices of locations (across the output field) to aggregate over. Defaults to NULL, which uses all outputs
#'
#' @return A plot of the emulator samples, with true replicates overlaid
#'
#' @export
PlotReplicates <- function(Samples, Replicates, locs = NULL){
  ell <- dim(Samples)[1]
  ns <- dim(Samples)[2]
  r <- ncol(Replicates)

  if (is.null(locs)){
    locs <- 1:ell
  }

  # Subsetting to chosen locations
  Samples <- Samples[locs,]
  Replicates <- Replicates[locs,]

  plot_data <- data.frame(Input = 1:length(locs),
                          Output = c(Samples),
                          s = rep(1:ns, each = length(locs)))
  plot_data$Output[plot_data$Output == 0] <- 0.1
  
  rep_data <- data.frame(Input = 1:length(locs),
                         Output = c(Replicates),
                         s = rep(ns + (1:r), each = length(locs)))
  rep_data$Output[rep_data$Output == 0] <- 0.1

  plot <- ggplot(plot_data, aes(.data$Input, .data$Output, col = as.factor(.data$s))) +
    geom_line() +
    geom_line(data = rep_data, aes(.data$Input, .data$Output, linetype = as.factor(.data$s)), size = 0.5, col = 'black') +
    scale_linetype_manual(values = rep(1, r)) +
    geom_line(data = data.frame(aggregate(Output ~ Input, plot_data, mean)), col = 'red', size = 1) +
    geom_line(data = data.frame(aggregate(Output ~ Input, plot_data, quantile, probs = 0.025)), col = 'red', size = 1, linetype = 'dashed') +
    geom_line(data = data.frame(aggregate(Output ~ Input, plot_data, quantile, probs = 0.975)), col = 'red', size = 1, linetype = 'dashed') +
    scale_colour_manual(values = rep('grey', max(plot_data$s))) +
    scale_y_log10(breaks = c(0.1,1,10,100,1000,10000,100000), labels = c('0','1','10','100','1000','10000','100000')) +
    theme_bw() +
    theme(legend.position = 'none')
  
  if (length(locs)<=15){
    plot <- plot + scale_x_continuous(breaks=1:length(locs))
  }

  return(plot)
}
