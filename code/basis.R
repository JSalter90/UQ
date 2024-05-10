# Functions relating to basis creation, manipulation

#' Formatting data
#'
#' Formats data so that it is in the correct form for use in other functions, and calculates the (weighted) SVD basis of the ensemble
#'
#' @param data a matrix containing individual fields in the columns (i.e. the matrix has dimension lxn)
#' @param weightinv the inverse of lxl positive definite weight matrix W. If NULL, the identity matrix is used
#' @param RemoveMean if TRUE, centres the data prior to calculating the basis
#' @param StoreEigen if TRUE, stores Q, lambda from eigendecomposition of W (in order to make later calculations more efficient)
#'
#' @return \item{tBasis}{The (weighted) SVD basis of the centred ensemble if RemoveMean = TRUE, of the original data otherwise}
#' \item{CentredField}{The centred data if RemoveMean = TRUE, the original data otherwise.}
#' \item{EnsembleMean}{The mean across the columns of the data. A zero vector if RemoveMean = FALSE}
#' \item{}
#'
#' @export
MakeDataBasis <- function(data, weightinv = NULL, W = NULL, RemoveMean = TRUE, StoreEigen = TRUE){
  if (RemoveMean == TRUE){
    EnsembleMean <- apply(data, 1, mean)
    CentredField <- 0*data
    for (i in 1:dim(data)[2]){
      CentredField[,i] <- data[,i] - EnsembleMean
    }
  }
  else {
    EnsembleMean <- c(rep(0, dim(data)[1]))
    CentredField <- data
  }
  if (is.null(W)){
    tSVD <- wsvd(t(CentredField), weightinv = weightinv)
    tBasis <- tSVD$v
    if (StoreEigen == TRUE){
      Q <- tSVD$Q
      Lambda <- tSVD$Lambda
      return(list(tBasis = tBasis, CentredField = CentredField, EnsembleMean = EnsembleMean, Q = Q, Lambda = Lambda))
    }
    else {
      return(list(tBasis = tBasis, CentredField = CentredField, EnsembleMean = EnsembleMean))
    }
  }
  else if (!is.null(W) & is.null(weightinv)){
    eig <- eigen(W)
    Q <- eig$vectors
    Lambda <- 1 / eig$values
    Winv <- Q %*% diag(Lambda) %*% t(Q)
    attr(Winv, 'diagonal') <- FALSE
    attr(Winv, 'identity') <- FALSE
    tSVD <- wsvd(t(CentredField), weightinv = Winv, Q = Q, Lambda = Lambda)
    tBasis <- tSVD$v
    if (StoreEigen == TRUE){
      return(list(tBasis = tBasis, CentredField = CentredField, EnsembleMean = EnsembleMean, Q = Q, Lambda = Lambda, Winv = Winv))
    }
    else {
      return(list(tBasis = tBasis, CentredField = CentredField, EnsembleMean = EnsembleMean, Winv = Winv))
    }
  }
}





#' Combine input data with basis projections
#' 
#' Processes data into the correct form for emulation of q coefficients
#' 
#' @param Design n x p matrix of inputs. Usually should already be scaled, and in the same order as the fields in DataBasis
#' @param DataBasis Object created by MakeDataBasis, containing basis, centred fields, and other information
#' @param q The number of basis vectors to project onto
#' @param Noise Whether to include a noise vector (sometimes used in selection of the GP mean function)
#' @param weightinv matrix to use for projection. By default, looks for this in DataBasis, can be overwritten in the function arguments if desire a different weighted projection
#' 
#' @return Inputs and outputs required for emulation
#' 
GetEmData <- function(Design, DataBasis, q = NULL, Noise = TRUE, weightinv = NULL){
  if(Noise){
    Noise <- runif(length(Design[,1]),-1,1)
    Design <- cbind(Design, Noise)
  }
  if (is.null(weightinv)){
    weightinv <- DataBasis$Winv
  }
  tData <- Project(DataBasis$CentredField, DataBasis$tBasis[,1:q], weightinv = weightinv)
  colnames(tData) <- paste0('C', 1:ncol(tData)) 
  tData <- data.frame(Design, tData)
  return(tData)
}

GetEmulatableData <- GetEmData


#' Weighted singular value decomposition
#'
#' Calculates the SVD basis across the output, given the inverse of W.
#'
#' @param data n x l matrix to calculate basis from (i.e. rows are output fields).
#' @param weightinv l x l inverse of W. If NULL, calculates standard SVD.
#' @param Q l x l matrix from eigen decomposition of W^{-1}, if already have this then speeds up calculation of basis
#' @param Lambda vector from eigen decomposition of W^{-1}, if already have this then speeds up calculation of basis
#'
#' @return The weighted SVD of the data.
#'
wsvd <- function(data, weightinv = NULL, Q = NULL, Lambda = NULL){
  if (is.null(weightinv)){
    svd_output <- svd(data)
  }
  else {
    stopifnot(dim(data)[2] == dim(weightinv)[1])
    if (is.null(Q) & attributes(weightinv)$diagonal == FALSE){
      eig <- eigen(weightinv)
      Q <- eig$vectors
      Lambda <- eig$values
      data_w <- data %*% Q %*% diag(sqrt(Lambda)) %*% t(Q)
      svd_output <- svd(data_w)
      svd_output$v <- t(t(svd_output$v) %*% Q %*% diag(1 / sqrt(Lambda)) %*% t(Q))
      svd_output$Q <- Q
      svd_output$Lambda <- Lambda
    }
    else if (is.null(Q) & attributes(weightinv)$diagonal == TRUE){
      diag_values <- diag(weightinv)
      data_w <- data %*% diag(sqrt(diag_values))
      svd_output <- svd(data_w)
      svd_output$v <- t(t(svd_output$v) %*% diag(1 / sqrt(diag_values)))
    }
    else if (!is.null(Q)){
      data_w <- data %*% Q %*% diag(sqrt(Lambda)) %*% t(Q)
      svd_output <- svd(data_w)
      svd_output$v <- t(t(svd_output$v) %*% Q %*% diag(1 / sqrt(Lambda)) %*% t(Q))
      svd_output$Q <- Q
      svd_output$Lambda <- Lambda
    }
  }
  return(svd_output)
}

#' Matrix inversion via cholesky decomposition
#'
#' Inverts matrix W, assigning attributes for whether W is diagonal, to speed up other calculations.
#'
#' @param W square positive definite variance matrix
#'
#' @return Inverse of W, with attributes 'identity' and 'diagonal', used by other functions in the package to make calculations more efficient.
#'
#' @examples Winv <- GetInverse(diag(100))
#' attributes(Winv) # diagonal = TRUE, identity = TRUE
#'
#' Winv2 <- GetInverse(runif(100,0.1,1)*diag(100))
#' attributes(Winv2) # diagonal = TRUE, identity = FALSE
#'
#' Winv3 <- GetInverse(seq(0.1,1,length=100) %*% t(seq(0.1,1,length=100)) + 0.1*diag(100))
#' attributes(Winv3) # diagonal = FALSE, identity = FALSE
#'
#' @export
GetInverse <- function(W){
  diagmat <- all(W[lower.tri(W)] == 0, W[upper.tri(W)] == 0)
  if (diagmat == TRUE){
    InvW <- diag(1 / diag(W))
  }
  else {
    Q <- chol(W)
    y <- backsolve(Q, diag(dim(W)[1]), transpose = TRUE)
    InvW <- crossprod(y, y)
  }
  attr(InvW, 'diagonal') <- diagmat
  if (all(diag(InvW) == 1) & diagmat == TRUE){
    attr(InvW, 'identity') <- TRUE
  }
  else {
    attr(InvW, 'identity') <- FALSE
  }
  return(InvW)
}


#' Projection onto a basis
#'
#' Calculates the coefficients given by projecting data onto a basis
#'
#' @param data Data matrix to be projected, where each column is a representation on the original field
#' @param basis Basis matrix
#' @param weightinv If NULL, uses standard SVD projection. Otherwise, uses weighted projection.
#'
#' @return Matrix of basis coefficients
#'
#' @examples # First generate some data
#'
#' l <- 100 # dimension of output
#' n <- 10 # number of runs
#' DataBasis <- MakeDataBasis(data = matrix(runif(l*n), nrow=l, ncol=n), RemoveMean = TRUE) # data is 100x10
#'
#' # Project the (centred) ensemble onto the first 3 vectors of the SVD basis
#'
#' Coefficients <- Project(data = DataBasis$CentredField, basis = DataBasis$tBasis[,1:3])
#'
#' # Instead of projecting using W = I, define a W with varying diagonal
#'
#' W <- runif(l, 1, 5) * diag(l) # 100x100 diagonal matrix
#' W_inv <- GetInverse(W) # inverse needed for projection
#' Coefficients_weighted <- Project(data = DataBasis$CentredField, basis = DataBasis$tBasis[,1:3], weightinv = W_inv)
#'
#' @export
Project <- function(data, basis, weightinv = NULL){
  d <- dim(data)[2]
  if (is.null(d)){
    d <- 1
  }
  p <- dim(basis)[2]
  l <- dim(basis)[1]
  if (is.null(p)){
    p <- 1
  }
  if (d == 1){
    data <- as.vector(data)
  }
  if (is.null(weightinv)){
    weightinv <- 0 # just need to set as something that isn't NULL so can give attribute
    attr(weightinv, 'diagonal') <- attr(weightinv, 'identity') <- TRUE
  }
  if (attributes(weightinv)$identity == TRUE){
    V <- t(basis) %*% basis
    Q <- chol(V)
    y <- backsolve(Q, diag(p), transpose = TRUE)
    x <- backsolve(Q, t(basis) %*% data, transpose = TRUE)
    scores <- crossprod(y, x)
  }
  else if (attributes(weightinv)$diagonal == TRUE) {
    V <- t(basis) %*% (diag(weightinv) * basis)
    Q <- chol(V)
    y <- backsolve(Q, diag(p), transpose = TRUE)
    tmp <- t(basis) %*% (diag(weightinv) * data)
    x <- backsolve(Q, tmp, transpose = TRUE)
    scores <- crossprod(y, x)
  }
  else {
    V <- t(basis) %*% weightinv %*% basis
    Q <- chol(V)
    y <- backsolve(Q, diag(p), transpose = TRUE)
    x <- backsolve(Q, t(basis) %*% weightinv %*% data, transpose = TRUE)
    scores <- crossprod(y, x)
  }
  # Assign coefficient labels to columns
  scores <- t(scores)
  colnames(scores)[1:p] <- paste("C",1:p,sep="")
  return(scores)
}

CalcScores <- Project



#' Field reconstructions from coefficients
#'
#' Given a vector of coefficients for a basis, calculates the field
#'
#' @param coeffs Coefficient vector
#' @param basis Basis matrix
#' @return Reconstructed field.
#'
#' @export
Recon <- function(coeffs, basis){
  if (is.null(dim(basis)[2])){
    q <- 1
  }
  else {
    q <- dim(basis)[2]
  }
  stopifnot(length(coeffs) == q)
  if (is.null(dim(basis)[2])){
    reconstruction <- basis*as.numeric(coeffs)
  }
  else {
    reconstruction <- basis%*%as.numeric(coeffs)
  }
  return(reconstruction)
}


#' Project and reconstruct a given field
#'
#' Gives the reconstruction of a field using a basis, by projecting and back-projecting on this basis.
#'
#' @param field Vector over original field
#' @param basis Basis matrix
#'
#' @return Reconstruction of the original field.
#'
#' @examples
#'
#' @export
ReconField <- function(field, basis, ...){
  nb <- is.null(dim(basis))
  if(!nb)
    basis1 <- basis[,1]
  else
    basis1 <- basis
  field <- c(field)
  mask <- which(is.na(field-basis1))
  if(length(mask)>0){
    recons <- rep(NA, length(field))
    field <- field[-mask]
    if(nb)
      basis <- basis[-mask]
    else
      basis <- basis[-mask,]
    proj <- Project(field, basis, ...)
    recons.partial <- Recon(proj, basis)
    recons[-mask] <- recons.partial
  }
  else{
    proj <- Project(field, basis, ...)
    recons <- Recon(proj, basis)
  }
  return(recons)
}

ReconObs <- ReconField


#' Matrix projection
#'
#' Projects a variance matrix onto a given basis
#'
#' @param mat A square matrix to be projected onto the basis
#' @param basis The basis to project with
#' @param weightinv The inverse of positive definite matrix W. If NULL, uses the standard projection, otherwise projects in the norm given by W.
#'
#' @return The projection of the original matrix on the basis.
#'
#' @export
VarProj <- function(mat, basis, weightinv = NULL){
  if (is.null(weightinv)){
    proj <- t(basis) %*% mat %*% basis
  }
  else {
    V <- t(basis) %*% weightinv %*% basis
    Q <- chol(V)
    y <- backsolve(Q, diag(dim(basis)[2]), transpose = TRUE)
    x <- backsolve(Q, t(basis) %*% weightinv, transpose = TRUE)
    comp <- crossprod(y, x)
    proj <- comp %*% mat %*% t(comp)
  }
  return(proj)
}


#' Reconstruction error
#'
#' Calculates the reconstruction error, R_W(basis, obs), of the observations given a basis and W.
#'
#' @param obs The observations
#' @param basis Basis to project and reconstruct the observations with
#' @param weightinv Inverse of weight matrix W. If NULL (default), calculates the mean squared error
#' @param scale If TRUE, scales by the dimension (so analogous to mean squared error)
#'
#' @return The reconstruction error
#'
#' @export
ReconError <- function(obs, basis, weightinv = NULL, scale = TRUE){
  if (is.null(weightinv)){
    weightinv <- 0
    attr(weightinv, 'diagonal') <- attr(weightinv, 'identity') <- TRUE
  }
  field <- ReconField(obs, basis, weightinv)
  A <- c(obs) - field
  mask <- which(is.na(A))
  if(length(mask)>0){
    A <- A[-mask]
  }
  if (scale == TRUE){
    s <- length(c(obs))-length(mask)
  }
  else {
    s <- 1
  }
  if (attributes(weightinv)$diagonal == FALSE){
    if(length(mask)>0){
      warning("Implicit assumption that weight specified on the full field even though applying a mask to missing obs/ensemble grid boxes")
      weightinv <- weightinv[-mask,-mask]
    }
    wmse <- (t(A) %*% weightinv %*% A)/ s
  }
  else {
    if (attributes(weightinv)$identity == TRUE){
      wmse <- crossprod(A)/ s
    }
    else {
      wmse <- crossprod(A/(1/diag(weightinv)), A)/ s
    }
  }
  return(as.numeric(wmse))
}




#' Reconstruction errors
#' 
#' Calculates reconstruction error of observations for each truncation of a basis
#' 
#' @param obs (Centred) observation vector
#' @param basis The basis to project onto
#' 
#' @return Vector of reconstruction errors for the observations
#' 
#' @export
errors <- function(basis, obs, weightinv=NULL){
  p <- dim(basis)[2]
  err <- numeric(p)
  if (is.null(weightinv)){
    weightinv <- diag(dim(basis)[1])
    attr(weightinv, 'diagonal') <- attr(weightinv, 'identity') <- TRUE
  }
  for (i in 1:p){
    err[i] <- ReconError(obs, basis[,1:i], weightinv)
  }
  return(err)
}


#' Calculating the proportion of data explained by a basis
#'
#' Calculates the proportion of the data that is explained by projection onto a basis.
#'
#' @param basis The basis
#' @param data The data to be explained
#' @param weightinv Inverse of W (identity if NULL)
#' @param total_sum The total sum of squares of the data with respect to W
#' @param psi t(original_basis) %*% weightinv %*% original_basis, where the new basis is a linear combination of some original basis
#' @param basis_lincom Vector of linear combinations (if new basis is a linear combination of some original basis)
#'
#' @return The proportion of variability in the data that is explained by the basis
#'
#' @export
VarExplained <- function(basis, data, weightinv = NULL, total_sum = NULL, psi = NULL, basis_lincom = NULL){
  coeffs <- t(Project(data, basis, weightinv))
  recon <- basis %*% coeffs
  if (is.null(weightinv)){
    explained <- crossprod(c(recon))/crossprod(c(data))
  }
  else {
    if (is.null(psi)){
      if (attributes(weightinv)$diagonal == TRUE){
        explained_num <- sum(t(recon)^2 %*% diag(weightinv))
      }
      else {
        explained_num <- sum(diag(t(recon) %*% weightinv %*% recon))
      }
    }
    else {
      stopifnot(!is.null(basis_lincom))
      explained_num <- t(coeffs) %*% t(basis_lincom) %*% psi %*%
        basis_lincom %*% coeffs
      explained_num <- sum(diag(explained_num))
    }
    #explained_num <- 0
    #for (i in 1:dim(data)[2]){
    #  explained_num <- explained_num + t(recon[,i]) %*% weightinv %*% recon[,i]
    #}
    #explained_den <- 0
    #for (i in 1:dim(data)[2]){
    #  explained_den <- explained_den + t(data[,i]) %*% weightinv %*% data[,i]
    #}
    if (is.null(total_sum)){
      if (attributes(weightinv)$diagonal == TRUE){
        explained_den <- sum(t(data)^2 %*% diag(weightinv))
      }
      else {
        explained_den <- sum(diag(t(data) %*% weightinv %*% data))
      }
    }
    else {
      explained_den <- total_sum
    }
    explained <- explained_num / explained_den
  }
  return(explained)
}


#' Number of basis vectors required to explain proportion of data
#'
#' Finds the truncated basis that explains a set proportion of the variability in the data.
#'
#' @param basis Basis matrix
#' @param data Data matrix
#' @param vtot The total proportion of variability in the data to be explained by the truncated basis
#' @param weightinv The inverse of W
#'
#' @return The number of basis vectors required to explain vtot of the data.
#'
#' @export
ExplainT <- function(DataBasis, vtot = 0.95, weightinv = NULL){
  v <- 0
  q <- 0
  while (v < vtot & q < dim(DataBasis$tBasis)[2]){
    v <- VarExplained(DataBasis$tBasis[,1:(q+1)], DataBasis$CentredField, weightinv)
    q <- q + 1
  }
  return(q)
}


#' Takes mean/variance from set of basis emulators and samples reconstructed fields
BasisEmSamples <- function(BasisPred, DataBasis, ns = 100, AddMean = TRUE, ReturnAll = TRUE, BasisUncertainty = TRUE, ...){
  n <- nrow(BasisPred$Expectation)
  q <- ncol(BasisPred$Expectation)
  if (is.null(n)){ # i.e. a single vector has been provided
    n <- 1
    q <- length(BasisPred$Expectation) 
  }
  
  ell <- dim(DataBasis$tBasis)[1]
  Basis <- DataBasis$tBasis[,1:q]
  
  if (AddMean){
    mu <- DataBasis$EnsembleMean
  }
  else {
    mu <- 0*DataBasis$EnsembleMean
  }
  
  if (BasisUncertainty){
    BasMinusQ <- DataBasis$tBasis[,-(1:q)] # get deleted vectors
    DeletedCoeffs <- Project(DataBasis$CentredField, BasMinusQ, ...) # project onto these vectors
    EstVar <- apply(DeletedCoeffs, 2, var) # variance on these vectors
    
    # Want full basis when reconstructing now
    Basis <- DataBasis$tBasis
    
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
    em_samp <- matrix(0, ell, ns)
    for (s in 1:ns){
      samp <- rnorm(q, 
                    mean = BasisPred$Expectation,
                    sd = sqrt(BasisPred$Variance))
      rec <- mu + Recon(samp, Basis)
      em_samp[,s] <- rec
    }
  }
  
  if (n > 1){
    em_samp <- array(0, dim = c(ell, ns, n))
    for (i in 1:n){
      for (s in 1:ns){
        samp <- rnorm(q, 
                      mean = BasisPred$Expectation[i,],
                      sd = sqrt(BasisPred$Variance[i,]))
        rec <- mu + Recon(samp, Basis)
        em_samp[,s,i] <- rec
      }
    }
  }
  
  if (ReturnAll){
    return(em_samp)
  }
  
  else {
    
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
    
    return(list(mean = samp_mean,
                lower = lower,
                upper = upper))
  }
}


#' This mostly does the same as MakeDataBasis, likely remove once added any additional functionality to MakeDataBasis
# Create basis (DataBasis object) and other preliminary calculations (inverting W)
# CreateBasis <- function(data, type = 'SVD', W = NULL, RemoveMean = TRUE){
#   # Invert W, creating Winv and tagging with attributes to enable fast calculations later on
#   if (is.null(W)){
#     W <- diag(nrow(data))
#   }
#   
#   Winv <- GetInverse(W)
#   
#   if (type %in% c('SVD', 'svd', 'L2')){
#     DataBasis <- MakeDataBasis(data, weightinv = NULL, RemoveMean = RemoveMean, StoreEigen = TRUE)
#   }
#   
#   if (type %in% c('WSVD', 'wsvd', 'wSVD')){
#     DataBasis <- MakeDataBasis(data, weightinv = Winv, RemoveMean = RemoveMean, StoreEigen = TRUE)
#   }
#   
#   DataBasis$Type <- type
#   
#   if (DataBasis$Type %in% c('svd', 'L2')){
#     DataBasis$Type <- 'SVD'
#   }
#   
#   if (DataBasis$Type %in% c('WSVD', 'wsvd')){
#     DataBasis$Type <- 'wSVD'
#   }
#   
#   return(DataBasis)
# }




#' Combines some other functions, review to see whether adds anything unique
# AssessBasis <- function(DataBasis, Obs){
#   
#   max_q <- dim(DataBasis$tBasis)[2]
#   PlotData <- data.frame(Size = 1:max_q, Error = numeric(max_q), Explained = numeric(max_q))
#   
#   if (!is.null(DataBasis$scaling)){
#     PlotData$Error <- errors(DataBasis$tBasis, Obs - DataBasis$EnsembleMean, DataBasis$Winv)*DataBasis$scaling^2
#   }
#   else {
#     PlotData$Error <- errors(DataBasis$tBasis, Obs - DataBasis$EnsembleMean, DataBasis$Winv)
#   }
#   
#   if (is.null(DataBasis$Winv)){
#     var_sum <- crossprod(c(DataBasis$CentredField))
#   }
#   else {
#     var_sum <- sum(diag(t(DataBasis$CentredField) %*% DataBasis$Winv %*% DataBasis$CentredField))
#   }
#   
#   for (i in 1:max_q){
#     PlotData$Explained[i] <- VarExplained(DataBasis$tBasis[,1:i], DataBasis$CentredField, DataBasis$Winv, total_sum = var_sum)
#   }
#   
#   PlotData$Explained <- round(PlotData$Explained, 10) # to make sure plots when = 1
#   
#   chi_bound <- qchisq(0.995, nrow(DataBasis$tBasis)) / nrow(DataBasis$tBasis)
#   max_y <- max(c(PlotData$Error, chi_bound + 0.05))
#   
#   var_plot <- ggplot(data = PlotData, aes(x = Size)) +
#     geom_line(aes(y = Error), col = 'red') +
#     geom_line(aes(y = Explained * max_y), col = 'blue') +
#     xlab('Basis size') +
#     scale_y_continuous(
#       #name = 'Error',
#       limits = c(0,max_y),
#       #breaks = seq(from = 0, to = 1*max_y, by = 0.25*max_y),
#       #labels = NULL,
#       sec.axis = sec_axis(trans=~./max_y, name="Explained")) +
#     #geom_vline(xintercept = q) +
#     geom_hline(yintercept = chi_bound, linetype = 'dashed', col = 'red') +
#     geom_hline(yintercept = 0.9*max_y, linetype = 'dashed', col = 'blue')
#   #theme(panel.grid.major = element_blank(), 
#   #      panel.grid.minor = element_blank())
#   
#   # Then don't need to do ExplainT, as already have the error/explained combinations
#   # Give list of (threshold, q)
#   thresholds <- c(0.8, 0.85, 0.9, 0.95, 0.99, 0.999)
#   q <- numeric(length(thresholds))
#   for (j in 1:length(q)){
#     q[j] <- which(PlotData$Explained >= thresholds[j])[1]
#   }
#   
#   return(list(plot = var_plot,
#               Errors = PlotData,
#               Truncations = data.frame(Threshold = thresholds,
#                                        q = q)))
# }



