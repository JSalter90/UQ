#' History matching
#'
#' Performing history matching, in 1D, for multiple independent outputs, or for high-dimensional fields, given observations, error (co)variances and corresponding emulators.
#'
#' This is a general function that can handle history matching across a wide variety of cases. Based on what is provided (some combination of emulators, sets of inputs to calculate the implausibility at, observations, ...), it recognises whether we are in the situation of multiple independent outputs or whether a single implausibility is to be calculated.
#' 
#' The possible cases, and the expected behaviour:
#' 
#' * `Obs` is a scalar, and so `emulator` should contain only a single emulator, and the 1D implausibility is calculated via `Impl1D`
#' * `Obs` is a vector of length k, and `emulator` is a list of k emulators. Each of the k outputs is treated independently, and the 1D implausibility is calculated via `Impl1D`
#' * `Obs` is a vector of length l, and `emulator` is a list of q < l emulators. Assumes that the leading q coefficients from `DataBasis$tBasis` have been emulated, and calculates `ImplField` (efficiently via basis decomposition)
#' 
#' Other cases (e.g., emulating l outputs independently, but wishing to calculate the l-dimensional implausiblity rather than l independent 1D implausibilities) are not supported by this function.
#'
#' @param emulator Either a single emulator, or a list of emulators, likely an object output by `FitEmulators`. If `NULL`,, then `TrueOutput` must be provided instead. 
#' @param design 
#' @param DataBasis If the multivariate implausibility is being calculated, need to provide the basis that the coefficents relate to. Likely an object output by `MakeDataBasis`.
#' @param Obs The observations, z.
#' @param ObsVar 
#' @param DiscVar 
#' @param TrueOutput 
#' @param PreviousWave 
#' @param ... 
#' 
#' @returns \item{Impl}{The implausibility (possible multiple) for each input, the size of NROY space, a logical vector, ...}
#' \item{Design}{The inputs that the implausibility was evaluated at, including a column indicating whether each input is in NROY.}
#' 
#' @export
HistoryMatch <- function(emulator = NULL, design = NULL, DataBasis = NULL, Obs, ObsVar, DiscVar = 0, TrueOutput = NULL, PreviousWave = NULL, ...){
  
  # If emulator not provided, must have been given TrueOutput, check
  if (is.null(emulator)){
    if (is.null(TrueOutput)){
      stop('Need to provide either an emulator or the true output')
    }
  }
  
  # In case a vector provided, convert to a matrix
  if (!(is.null(TrueOutput))){
    TrueOutput <- as.matrix(TrueOutput)
  }

  # If only 1 observation, check only been given a single emulator (or single true output)
  if (length(Obs) == 1){
    if (!(is.null(emulator))){
      if (is.null(emulator$method)){
        stop('1 observation provided, but either too many emulators, or emulator in incorrect format')
      }
    }
    else {
      if (ncol(TrueOutput) > 1){
        stop('1 observation provided, but more than 1 outputs')
      }
    }
    Multivariate <- FALSE
  }
  # If more than 1 obs, we might have only a single emulator (single coefficient)
  
  # Determine whether HM in 1D across multiple metrics, or have basis emulator
  else {
    if (!(is.null(emulator))){
      Multivariate <- ifelse(length(emulator) == length(Obs), FALSE, TRUE)
    }
    
    else {
      Multivariate <- ifelse(ncol(TrueOutput) == length(Obs), FALSE, TRUE)
    }
  }
  
  if (Multivariate & is.null(DataBasis)){
    stop('Mismatch between number of observations and either number of emulators or number of true outputs. 
         Either ensure these match up (1D HM), or provide a DataBasis object to do multivariate HM')
  }
  
  if (is.null(design) & is.null(PreviousWave)){
    stop('Requires design if there is no PreviousWave to inherit this from')
  }
  
  # Find wave number (for clear labelling of NROY spaces when multiple waves)
  if (is.null(PreviousWave)){
    wave <- 1
  }
  
  else {
    wave <- PreviousWave$wave + 1
  }

  # If previous HM wave provided, select points currently in NROY, and only predict for these
  if (!(is.null(PreviousWave))){
    current_nroy <- paste0('NROY',wave-1)
    nroy_inds <- which(PreviousWave$Design[,current_nroy])
    design <- PreviousWave$Design[nroy_inds,]
  }

  # Evaluate emulators at design
  if (is.null(TrueOutput)){
    preds <- Predict(emulator, design)
  }
  
  # If TrueOutput was provided, instead calculate implausibility based on truth (i.e., zero variance)
  else {
    preds <- NULL
    preds$Expectation <- TrueOutput
    preds$Variance <- 0*TrueOutput
  }

  if (Multivariate){
    impl <- ImplField(preds, DataBasis, Obs = Obs, ObsVar = ObsVar, DiscVar = DiscVar, ...)
  }
  
  else {
    impl <- Impl(preds, Obs = Obs, ObsVar = ObsVar, DiscVar = DiscVar, ...)
  }
  
  design <- as.data.frame(design)
  new_nroy <- paste0('NROY', wave)
  design[,new_nroy] <- impl$inNROY
  
  # Combine new NROY information with previous if it exists
  if (!(is.null(PreviousWave))){
    PreviousWave$Design[,new_nroy] <- FALSE
    PreviousWave$Design[nroy_inds[impl$inNROY],new_nroy] <- TRUE
    design <- PreviousWave$Design
  }

  # If there was a previous wave, update $NROY to give accurate percentage of original space
  if (!(is.null(PreviousWave))){
    impl$NROY <- impl$NROY * PreviousWave$NROY
  }
  
  if (Multivariate){
    return(list(Impl = impl$Impl, bound = impl$bound, NROY = impl$NROY, inNROY = impl$inNROY,
                Design = design,
                Preds = preds,
                wave = wave))
  }
  
  else {
    return(list(Impl = impl$Impl, NROY = impl$NROY, Pass_kth = impl$Pass_kth, Pass_N = impl$Pass_N, inNROY = impl$inNROY,
                Design = design,
                Preds = preds,
                wave = wave))
  }
}


#' High-dimensional implausibility
#'
#' Calculates the implausibility (for the \ell-dimensional field) for a sample of expectations and variances from basis emulators
#'
#' @param DataBasis object containing the basis used in emulation ($tBasis)
#' @param Obs The observations z, a vector of length l. Should not be centered as this happens internally if required, according to `DataBasis$EnsembleMean`.
#' @param Predictions Containing components `$Expectation` and `$Variance`, a matrix containing emulator expectations (variances), where a given row contains the expectations for the q emulated basis vectors, for some x
#' @param ObsVar observation error variance matrix
#' @param DiscVar discrepancy variance matrix. Defaults to NULL (zero), as can be included as part of `Error`
#' @param obs_inds If the observation vector is a subset of the full output, selects required indices from full output that correspond to observations. Defaults to NULL (all outputs observed)
#' @param weightinv if not NULL, the inverse of W = var_err + var_disc, used for projection
#' #' 
#' @returns \item{impl}{Vector of implausibilities corresponding to the rows of Expectation and Variance}
#' \item{bound}{The chi-squared bound for an \ell-dimensional field}
#' \item{nroy}{Proportion of parameter settings that are not ruled out, using bound}
#' \item{inNROY}{Vector indicating whether a parameter setting is ruled out}
#'
#' @export
ImplField <- function(Predictions, DataBasis, Obs, ObsVar, DiscVar = 0, obs_inds = NULL, weightinv = NULL, BasisUncertainty = TRUE){
  
  Expectation <- Predictions$Expectation
  Variance <- Predictions$Variance
  
  n <- nrow(Expectation) # number of prediction points
  q <- ncol(Expectation) # number of emulated vectors
  Basis <- DataBasis$tBasis[,1:q]
  ell <- nrow(Basis)

  if (!(q == ncol(Variance))){
    stop('Incorrect number of columns in Variance')
  }
  
  # If a subset is not selected, use all outputs
  if (is.null(obs_inds)){
    obs_inds <- 1:ell
  }

  # Check provided observation vector is consistent with the length of the model subset selected 
  if (!(length(Obs) == length(obs_inds))){
    stop('The observation vector is not the same length as the model output')
  }
  
  # Centre the observations (if DataBasis was centred)
  Obs <- Obs - DataBasis$EnsembleMean[obs_inds]
  
  # If provided single variance(s), assume iid
  if (length(ObsVar) == 1){
    ObsVar <- ObsVar*diag(length(obs_inds))
  }

  if (length(DiscVar) == 1){
    DiscVar <- DiscVar*diag(length(obs_inds))
  }
  
  W <- ObsVar + DiscVar
  if (is.null(weightinv)){
    weightinv <- GetInverse(W)
  }

  # Add uncertainty from discarded basis vectors?
  if (BasisUncertainty == TRUE){
    # These need to be projected in the same way that the emulated coefficients were
    if (is.null(DataBasis$Winv)){
      BasisVar <- DiscardedBasisVariance(DataBasis, q)
    }
    else {
      BasisVar <- DiscardedBasisVariance(DataBasis, q, weightinv = DataBasis$Winv)
    }
    W <- W + BasisVar[obs_inds,obs_inds] # basis variance exists over full output, but if we've got partial obs, only need these rows/columns
    weightinv <- GetInverse(W) # need to re-define W^-1 to include this extra variance, to enable fast calculation of full I(x)
  }
  
  # A faster implementation, if the full basis is still relatively small, and we are not considering a large number of points
  # if (BasisUncertainty == TRUE){
  #   # These need to be projected in the same way that the emulated coefficients were
  #   BasMinusQ <- DataBasis$tBasis[,-(1:q)]
  #   DeletedCoeffs <- Project(DataBasis$CentredField, BasMinusQ)
  #   EstVar <- apply(DeletedCoeffs, 2, var) # vector of variances for deleted vectors
  #   
  #   # Re-define basis as full basis
  #   Basis <- DataBasis$tBasis
  #   
  #   # Append columns of zeros to expectation to represent unemulated directions
  #   Expectation <- cbind(Expectation, matrix(0, ns, ncol(Basis) - q))
  #   
  #   # Append columns with variance for each deleted direction
  #   Variance <- cbind(Variance, matrix(rep(EstVar, each = ns), ns, ncol(Basis) - q))
  # }
  
  R_W <- ReconError(Obs, Basis[obs_inds,], weightinv = weightinv)
  
  # Project observations onto basis if required
  if (length(Obs) == length(obs_inds)){
    ObsProj <- Project(Obs, Basis[obs_inds,], weightinv = weightinv)
  }
  
  # Project variance matrices onto basis if required
  if (nrow(ObsVar) == length(obs_inds)){
    WProj <- ProjectVar(W, Basis[obs_inds,], weightinv = weightinv)
  }
  impl <- as.numeric(parallel::mclapply(1:n, 
                                        function(i) ImplCoeff(Expectation[i,], Variance[i,], ObsProj, WProj, 0*WProj),
                                        mc.cores = n_cores))
  impl <- impl + rep(R_W, n) 
  bound <- stats::qchisq(0.995, length(obs_inds))
  nroy <- sum(impl < bound)/n
  inNROY <- impl < bound

  return(list(Impl = impl, 
              bound = bound, 
              NROY = nroy, 
              inNROY = inNROY))
}





#' Coefficient implausibility
#'
#' Calculates the coefficient implausibility for a single x, given projected quantities
#'
#' @param Expectation length q vector with emulator expectations.
#' @param Variance length q vector with emulator variances.
#' @param Obs vector with projected observations, must be the same length as `Expectation` and `Variance`.
#' @param ObsVar projected observation error variance matrix, must have dimension qxq.
#' @param DiscVar projected discrepancy variance matrix, must have dimension qxq.
#'  
#' @return The coefficient implausibility (given the matrix used in projection)
#'
#' @export
ImplCoeff <- function(Expectation, Variance, Obs, ObsVar, DiscVar){
  
  if (!(length(Expectation) == length(Variance))){
    stop('Emulator expectation and variance have different dimensions')
  }
  
  if (!(length(Expectation) == length(Obs))){
    stop('Emulator expectation and projected observations have different dimensions')
  }
  
  if (!(length(Expectation) == nrow(ObsVar))){
    stop('Emulator expectation and projected observation error matrix have inconsistent dimensions')
  }
  
  if (!(length(Expectation) == nrow(DiscVar))){
    stop('Emulator expectation and projected discrepancy matrix have inconsistent dimensions')
  }
  
  V <- ObsVar + DiscVar + diag(Variance)
  Q <- chol(V)
  proj.output <- Expectation
  y <- backsolve(Q, as.vector(Obs - proj.output), transpose = TRUE)
  impl <- crossprod(y,y)
  return(impl)
}




#' High-dimensional implausibility without emulation
#'
#' @param Fields 
#' @param Obs 
#' @param ObsVar 
#' @param DiscVar 
#'
#' @return
#' @export
#'
#' @examples
ImplFieldTrue <- function(Fields, Obs, ObsVar, DiscVar = 0){

  ell <- nrow(Fields)
  n <- ncol(Fields)
  
  if (!(ell == length(Obs))){
    stop('Observations and Fields have different lengths')
  }
  
  # If provided single variance(s), assume iid
  if (length(ObsVar) == 1){
    ObsVar <- ObsVar*diag(ell)
  }
  
  if (length(DiscVar) == 1){
    DiscVar <- DiscVar*diag(ell)
  }
  
  tmp_inverse <- GetInverse(ObsVar + DiscVar)
  
  Impl <- rep(NA, n)
  for (i in 1:n){
    Impl[i] <- t(Obs - Fields[,i]) %*% tmp_inverse %*% c(Obs - Fields[,i])
  }
  
  return(Impl)
}






#' Instead of calculating via coeff + recon_error, could reconstruct everything
#' Should be identical
#' For testing only
#'
#' @param Predictions 
#' @param DataBasis 
#' @param Obs 
#' @param ObsVar 
#' @param DiscVar 
ImplFieldRecon <- function(Predictions, DataBasis, Obs, ObsVar, DiscVar = 0){
  
  Expectation <- Predictions$Expectation
  Variance <- Predictions$Variance
  
  n <- nrow(Expectation)
  q <- ncol(Expectation)
  ell <- nrow(DataBasis$tBasis)
  
  if (!(q == ncol(Variance))){
    stop('Incorrect number of columns in Variance')
  }
  
  # If provided single variance(s), assume iid
  if (length(ObsVar) == 1){
    ObsVar <- ObsVar*diag(ell)
  }
  
  if (length(DiscVar) == 1){
    DiscVar <- DiscVar*diag(ell)
  }
  
  BasisVar <- DiscardedBasisVariance(DataBasis, q)
  
  Impl <- rep(NA, n)
  for (i in 1:n){
    recon_mean <- Recon(Expectation[i,], DataBasis$tBasis[,1:q]) + DataBasis$EnsembleMean
    recon_var <- DataBasis$tBasis[,1:q] %*% diag(Variance[i,]) %*% t(DataBasis$tBasis[,1:q])
    Impl[i] <- t(Obs - recon_mean) %*% solve(recon_var + ObsVar + DiscVar + BasisVar) %*% c(Obs - recon_mean)
  }
  
  return(Impl)
}






#' Setting discrepancy multiple
#'
#' Scaling the discrepancy to ensure that the observations won't be ruled out
#'
#' @param Basis full basis
#' @param q where the basis is truncated
#' @param obs observations
#' @param level quantile of the chi-squared distribution to use (< 0.995)
#' @param weightinv inverse of W, to use in projection
#'  
#' @return scalar to be used as discrepancy multiplier, to ensure observations not ruled out
#'
#' @export
SetDiscrepancy <- function(Basis, q, obs, level = 0.95, weightinv = NULL){
  TruncatedError <- ReconError(obs, Basis[,1:q], weightinv = weightinv)
  ell <- dim(Basis)[1]
  b <- stats::qchisq(level, ell)
  DiscMultiplier <- c(TruncatedError / b)
  return(DiscMultiplier)
}


#' Adding discarded basis variance
#' 
#' Function that gives variance given by deleted basis vectors (Wilkinson 2010), used within HistoryMatch
#' 
#' @param DataBasis - object containing basis, centred field etc.
#' @param q - where the basis is truncated
#' @param weightinv inverse of W, to use in projection
#' 
#' @return Matrix containing uncertainty due to basis truncation
#' 
#' @export
DiscardedBasisVariance <- function(DataBasis, q, weightinv = NULL){
  BasMinusQ <- DataBasis$tBasis[,-(1:q)]
  DeletedCoeffs <- Project(DataBasis$CentredField, BasMinusQ, weightinv)
  EstVar <- apply(DeletedCoeffs, 2, stats::var)
  DeletedBasisVar <- BasMinusQ %*% diag(EstVar) %*% t(BasMinusQ)
  return(DeletedBasisVar)
}





#' Title
#'
#' @param x 
#' @param k 
#'
#' @returns
#' @export
#'
#' @examples
kth_max <- function(x,k) {
  sorted_values <- sort(x, decreasing = TRUE)
  sorted_values[k]  # returns the kth maximum
}




# How many pass kth constraint
#' Title
#'
#' @param impl 
#' @param bound 
#'
#' @returns
#' @export
#'
#' @examples
pass_kth_metric <- function(impl, bound = 3){
  ns <- nrow(impl)
  apply(impl < bound, 2, sum) / ns
}




# How many constraints each run passes. RENAME? pass_by_run?
#' Title
#'
#' @param impl 
#' @param bound 
#'
#' @returns
#' @export
#'
#' @examples
pass_metric <- function(impl, bound = 3){
  ns <- nrow(impl)
  summary(as.factor(apply(impl < bound, 1, sum, na.rm = TRUE))) / ns
}



#' Implausibility
#'
#' 1D impl for 1 or more metrics
#'
#' @param em_pred emulator predictions (containing `$Expectation`, `$Variance`)
#' @param obs_error vector of observation error variances corresponding to each output
#' @param kmax kth max implausibility, defaults to 1 (i.e. take max implausibility across outputs)
#'
#' @return \item{overall_impl}{Overall implausibility for each of the left-out runs}
#' \item{}
#'
#' @export
Impl <- function(Predictions, Obs, ObsVar, DiscVar = 0, kmax = 1, bound = 3){

  Expectation <- Predictions$Expectation
  Variance <- Predictions$Variance
  
  # Number of outputs
  ell <- ncol(Expectation)
  
  # Number of prediction points
  n <- nrow(Expectation)

  if (is.null(ell)){
    ell <- 1
  }
  
  if (!(length(Obs) == ell)){
    stop('Inconsistency between number of emulator predictions and number of observations provided')
  }
  
  if (!(length(ObsVar) %in% c(1,ell))){
    stop('Inconsistency between number of emulator predictions and number of observation error variances provided')
  }
  
  if (length(ObsVar) == 1 & !(ell == 1)){
    ObsVar <- rep(ObsVar, ell)
  }
  
  if (!(length(DiscVar) %in% c(1,ell))){
    stop('Inconsistency between number of emulator predictions and number of discrepancy variances provided')
  }
  
  if (length(DiscVar) == 1 & !(ell == 1)){
    DiscVar <- rep(DiscVar, ell)
  }
  
  # Converting to vector
  Obs <- as.numeric(Obs)
  
  # Calculate implausibility
  impl <- parallel::mclapply(1:ell, 
                             function(k) Impl1D(Expectation[,k], Variance[,k], Obs[k], ObsVar[k], DiscVar[k]),
                             mc.cores = n_cores)
  impl <- matrix(unlist(impl), n, ell)

  # Proportion of samples that pass each metric individually
  pass_kth <- pass_kth_metric(impl, bound)
  
  # Number of constraints passed by each sample
  pass_N <- pass_metric(impl, bound)
  
  # Defining NROY based on kmax
  impl_kmax <- apply(impl, 1, kth_max, k = kmax)
  nroy <- mean(impl_kmax < bound)
  inNROY <- impl_kmax < bound

  return(list(Impl = impl,
              NROY = nroy,
              Pass_kth = pass_kth,
              Pass_N = pass_N,
              inNROY = inNROY))
}






#' Title
#'
#' @param expectation 
#' @param variance 
#' @param Obs 
#' @param ObsVar 
#' @param DiscVar 
#'
#' @returns
#' @export
#'
#' @examples
Impl1D <- function(expectation, variance, Obs, ObsVar, DiscVar = 0){
  abs(expectation - Obs) / sqrt(variance + ObsVar + DiscVar)
}






#' Title
#'
#' @param CandidatePoints 
#' @param n_new 
#' @param reps 
#'
#' @returns
#' @export
#'
#' @examples
NewDesign <- function(CandidatePoints, n_new, reps = 100){
  
  # If final column(s) of CandidatePoints are NROY indicator, remove
  nroy_inds <- grep('NROY', colnames(CandidatePoints))
  
  if (length(nroy_inds) > 0){
    CandidatePoints <- CandidatePoints[,-nroy_inds]
  }
  
  # Generate candidate designs of size n_new
  inds_list <- NULL
  dist_min <- numeric(reps)
  for (i in 1:reps){
    inds_list[[i]] <- sample(1:nrow(CandidatePoints), n_new)
    tmp_design <- CandidatePoints[inds_list[[i]],]
    dists <- fields::rdist(tmp_design)
    diag(dists) <- NA
    dist_min[i] <- min(c(dists), na.rm = TRUE)
    rm(dists)
  }
  
  k <- which.max(dist_min)
  
  new_design <- CandidatePoints[inds_list[[k]],]
  
  rownames(new_design) <- NULL
  
  return(list(design = new_design,
              inds = inds_list[[k]]))
}






#### not edited ####

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



#' Pairs plot of NROY space
#'
#' @param Design Data frame of inputs, possibly with a TRUE/FALSE final column named `NROY`.
#' @param k Indices relating to which columns of `Design` to plot. If NULL, plots all.
#' @param NROY Vector of TRUE/FALSE labels, corresponding to classification of `Design`. If NULL, uses the final column of `Design`.
#' @param Truth 
#' @param size Controlling size of points on lower half of plot
#'
#' @return
#' @export
#'
#' @examples
PlotNROY <- function(Design, k = NULL, NROY = NULL, Truth = NULL){
  # If a vector of NROY labels is not provided, assume this is the final column of Design
  if (is.null(NROY)){
    colnames(Design)[ncol(Design)] <- 'NROY'
  }
  
  else {
    Design <- data.frame(Design, NROY = NROY)
  }
  
  if (is.null(k)){
    k <- 1:(ncol(Design)-1)
  }
  
  size <- ifelse(nrow(Design) < 1000, 2, 0.5)
  
  # If only 2 inputs provided, single pairs plot
  if (length(k) == 2){
    Design$input1 <- Design[,k[1]]
    Design$input2 <- Design[,k[2]]
    
    plot <- ggplot(Design, aes(.data$input1, .data$input2, col = .data$NROY)) + 
      geom_point(size = size) +
      scale_color_manual(values = cols_validate[-1]) +
      labs(x = colnames(Design)[k[1]], y = colnames(Design)[k[2]]) +
      theme_bw()
    
    if (!is.null(Truth)){
      Truth <- as.data.frame(Truth)
      Truth$input1 <- Truth[,k[1]]
      Truth$input2 <- Truth[,k[2]]
      plot <- plot + geom_point(data = Truth, col = 'red', shape = 17, size = 4)
    }
  }
  
  else {
    plot <- GGally::ggpairs(Design, columns = k, aes(col = .data$NROY), 
                            upper = list(continuous = GGally::wrap("density", alpha = 0.5), combo = "box_no_facet"),
                            lower = list(continuous = GGally::wrap("points", alpha = 0.3, size = size), combo = GGally::wrap("dot_no_facet", alpha = 0.4)),
                            diag = list(continuous = GGally::wrap("densityDiag", alpha = 0.3)), 
                            legend = 1) +
      theme_bw() +
      theme(legend.position = "bottom") + 
      scale_colour_manual(values = cols_validate[-1]) +
      scale_fill_manual(values = cols_validate[-1])
    
    # Add truth
    if (!is.null(Truth)){
      Truth$NROY <- NA
      for(i in 2:length(k)){
        for (j in 1:(i-1)){
          tmp <- Truth
          tmp$input1 <- Truth[,k[j]]
          tmp$input2 <- Truth[,k[i]]
          p1 <- GGally::getPlot(plot, i, j) + geom_point(data = tmp, aes(x = .data$input1, y = .data$input2), col = 'red', shape = 17, size = 4)
          plot <- GGally::putPlot(plot, p1, i, j)
        }
      }
    }
  }
  return(plot)
}









#### NOT EDITED ####

# PlotEmulator <- function(impl, design, variance = FALSE){
#   thetapred <- unit_to_theta(impl[,c('X1','X2')])
#   impl$lambda <- thetapred$lambda
#   impl$C <- thetapred$C
#   
#   if (!(is.null(dim(design)))){
#     design <- unit_to_theta(design)
#   }
#   
#   if (variance){
#     plot <- ggplot(impl, aes(lambda, C, col = Var)) + 
#       geom_point() +
#       geom_point(data = data.frame(C = design$C, lambda = design$lambda), col = 'black', size = 5, shape = 15) +
#       scale_colour_viridis_c(option = 'H') +
#       theme_minimal(base_size = 14) +
#       theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"))
#   }
#   
#   else{
#     plot <- ggplot(impl, aes(lambda, C, col = Mean)) + 
#       geom_point() +
#       geom_point(data = data.frame(C = design$C, lambda = design$lambda), col = 'black', size = 5, shape = 15) +
#       scale_colour_viridis_c() +
#       theme_minimal(base_size = 14) +
#       theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"))
#   }
#   
#   return(plot)
# }
