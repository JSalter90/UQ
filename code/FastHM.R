#source('rotation_functions.R')

#' History matching
#'
#' Calculates the implausibility (for the \ell-dimensional field) for a sample of expectations and variances from basis emulators
#'
#' @param DataBasis object containing the basis used in emulation ($tBasis)
#' @param Obs observation vector (length \ell), must be centred
#' @param Expectation a matrix containing emulator expectations, where a given row contains the expectations for the q emulated basis vectors, for some x
#' @param Variance a matrix containing emulator variances, where a given row contains the variances for the q emulated basis vectors, for some x
#' @param Error observation error variance matrix
#' @param Disc discrepancy variance matrix
#' @param weightinv if not NULL, the inverse of W = var_err + var_disc, used for projection
#' #' 
#' @return \item{impl}{Vector of implausibilities corresponding to the rows of Expectation and Variance}
#' \item{bound}{The chi-squared bound for an \ell-dimensional field}
#' \item{nroy}{Proportion of parameter settings that are not ruled out, using bound}
#' \item{inNROY}{Vector indicating whether a parameter setting is ruled out}
#'
#' @export
HistoryMatch <- function(DataBasis, Obs, Expectation, Variance, Error, Disc, weightinv = NULL, BasisUncertainty = TRUE){
  q <- dim(Expectation)[2]
  Basis <- DataBasis$tBasis[,1:q]
  l <- dim(Basis)[1]
  stopifnot(q == dim(Variance)[2])
  W <- Error + Disc
  if (is.null(weightinv)){
    weightinv <- GetInverse(W)
  }
  nn <- dim(Expectation)[1]
  # Add uncertainty from discarded basis vectors?
  if (BasisUncertainty == TRUE){
    # These need to be projected in the same way that the emulated coefficients were
    if (is.null(DataBasis$Winv)){
      BasisVar <- DiscardedBasisVariance(DataBasis, q)
    }
    else {
      BasisVar <- DiscardedBasisVariance(DataBasis, q, weightinv = DataBasis$Winv)
    }
    W <- W + BasisVar
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
  #   Expectation <- cbind(Expectation, matrix(0, nn, ncol(Basis) - q))
  #   
  #   # Append columns with variance for each deleted direction
  #   Variance <- cbind(Variance, matrix(rep(EstVar, each = nn), nn, ncol(Basis) - q))
  # }

  R_W <- ReconError(Obs, Basis, weightinv = weightinv, scale = FALSE)
  # Project observations onto basis if required
  if (length(Obs) == l){
    ObsProj <- Project(Obs, Basis, weightinv = weightinv)
  }
  # Project variance matrices onto basis if required
  if (dim(Disc)[1] == l){
    WProj <- VarProj(W, Basis, weightinv = weightinv)
  }
  impl <- as.numeric(parallel::mclapply(1:nn, function(i) ImplCoeff(Expectation[i,], Variance[i,], ObsProj, WProj, 0*WProj)))
  impl <- impl + rep(R_W, nn) 
  bound <- qchisq(0.995, l)
  nroy <- sum(impl < bound)/nn
  inNROY <- impl < bound
  return(list(impl = impl, bound = bound, nroy = nroy, inNROY = inNROY))
}

#' Coefficient implausibility
#'
#' Calculates the coefficient implausibility for a single x, given projected quantities
#'
#' @param Expectation length q vector with emulator expectations
#' @param Variance length q vector with emulator variances
#' @param Obs projected observations
#' @param Error projected observation error variance matrix
#' @param Disc projected discrepancy variance matrix
#'  
#' @return The coefficient implausibility (given the matrix used in projection)
#'
#' @export
ImplCoeff <- function(Expectation, Variance, Obs, Error, Disc){
  V <- Error + Disc + diag(Variance)
  Q <- chol(V)
  proj.output <- Expectation
  y <- backsolve(Q, as.vector(Obs - proj.output), transpose = TRUE)
  impl <- crossprod(y,y)
  return(impl)
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
SetDiscrepancy <- function(tBasis, q, obs, level = 0.95, weightinv = NULL){
  TruncatedError <- ReconError(obs, tBasis[,1:q], weightinv = weightinv, scale = FALSE)
  l <- dim(tBasis)[1]
  b <- qchisq(level, l)
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
  EstVar <- apply(DeletedCoeffs, 2, var)
  DeletedBasisVar <- BasMinusQ %*% diag(EstVar) %*% t(BasMinusQ)
  return(DeletedBasisVar)
}
