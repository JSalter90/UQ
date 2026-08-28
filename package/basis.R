#' Functions for basis construction, manipulation

# Edit so that the input values/name are already included in the basis

#' Formatting data
#'
#' Formats data so that it is in the correct form for use in other functions, and calculates the (weighted) SVD basis of the ensemble
#'
#' @param data a matrix containing individual fields in the columns (i.e. the matrix has dimension lxn)
#' @param weightinv the inverse of lxl positive definite weight matrix W. If NULL, the identity matrix is used
#' @param RemoveMean if TRUE, centres the data prior to calculating the basis
#' @param StoreEigen if TRUE, stores Q, lambda from eigendecomposition of W (in order to make later calculations more efficient)
#'
#' @returns \item{tBasis}{The (weighted) SVD basis of the centred ensemble if RemoveMean = TRUE, of the original data otherwise}
#' \item{CentredField}{The centred data if RemoveMean = TRUE, the original data otherwise.}
#' \item{EnsembleMean}{The mean across the columns of the data. A zero vector if RemoveMean = FALSE}
#' \item{}
#'
#' @export
MakeDataBasis <- function(data, weightinv = NULL, W = NULL, RemoveMean = TRUE, Coeffs = TRUE, StoreEigen = TRUE){
  if (RemoveMean == TRUE){
    EnsembleMean <- apply(data, 1, mean)
    CentredField <- 0*data
    for (i in 1:ncol(data)){
      CentredField[,i] <- data[,i] - EnsembleMean
    }
  }
  else {
    EnsembleMean <- c(rep(0, nrow(data)))
    CentredField <- data
  }
  if (is.null(W)){
    tSVD <- wsvd(t(CentredField), weightinv = weightinv)
    tBasis <- tSVD$v
    
    # Remove final vector
    if (RemoveMean){
      tBasis <- tBasis[,-ncol(tBasis)]
    }
    
    tCoeffs <- NULL
    if (Coeffs){
      tCoeffs <- Project(CentredField, tBasis, weightinv = weightinv)
      colnames(tCoeffs) <- paste0('C', 1:ncol(tCoeffs))
    }

    if (StoreEigen & !(is.null(tSVD$Q))){
      Q <- tSVD$Q
      Lambda <- tSVD$Lambda
      return(list(tBasis = tBasis, CentredField = CentredField, EnsembleMean = EnsembleMean, Coeffs = tCoeffs, Q = Q, Lambda = Lambda))
    }
    else {
      return(list(tBasis = tBasis, CentredField = CentredField, EnsembleMean = EnsembleMean, Coeffs = tCoeffs))
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
    
    # Remove final vector
    if (RemoveMean){
      tBasis <- tBasis[,-ncol(tBasis)]
    }
    
    # Calculate projection
    tCoeffs <- NULL
    if (Coeffs){
      tCoeffs <- Project(CentredField, tBasis, weightinv = Winv)
      colnames(tCoeffs) <- paste0('C', 1:ncol(tCoeffs))
    }
    
    if (StoreEigen){
      return(list(tBasis = tBasis, CentredField = CentredField, EnsembleMean = EnsembleMean, Coeffs = tCoeffs, Q = Q, Lambda = Lambda, Winv = Winv))
    }
    else {
      return(list(tBasis = tBasis, CentredField = CentredField, EnsembleMean = EnsembleMean, Winv = Winv, Coeffs = tCoeffs))
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
#' @returns Inputs and outputs required for emulation
#' 
GetEmData <- function(Design, DataBasis, q = NULL, Noise = FALSE, weightinv = NULL, input_range = c(-1,1)){
  if (q > ncol(DataBasis$tBasis)){
    stop('q is higher than the number of basis vectors, reduce')
  }

  if (Noise){
    Noise <- stats::runif(nrow(Design),input_range[1],input_range[2])
    Design <- cbind(Design, Noise)
  }
  
  if (is.null(weightinv)){
    weightinv <- DataBasis$Winv
  }
  
  if (is.null(DataBasis$Coeffs)){
    tData <- Project(DataBasis$CentredField, DataBasis$tBasis[,1:q], weightinv = weightinv)
    colnames(tData) <- paste0('C', 1:ncol(tData)) 
  }
  
  else {
    tData <- DataBasis$Coeffs[,1:q]
  }
  
  # Combine inputs with projections
  tData <- data.frame(Design, tData)
  return(tData)
}




#' Weighted singular value decomposition
#'
#' Calculates the SVD basis across the output, given the inverse of W.
#'
#' @param data n x l matrix to calculate basis from (i.e. rows are output fields).
#' @param weightinv l x l inverse of W. If NULL, calculates standard SVD.
#' @param Q l x l matrix from eigen decomposition of W^{-1}, if already have this then speeds up calculation of basis
#' @param Lambda vector from eigen decomposition of W^{-1}, if already have this then speeds up calculation of basis
#'
#' @returns The weighted SVD of the data.
#'
wsvd <- function(data, weightinv = NULL, Q = NULL, Lambda = NULL){
  if (is.null(weightinv)){
    svd_output <- svd(data)
  }
  else {
    stopifnot(ncol(data) == nrow(weightinv))
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
#' @returns Inverse of W, with attributes 'identity' and 'diagonal', used by other functions in the package to make calculations more efficient.
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
    y <- backsolve(Q, diag(nrow(W)), transpose = TRUE)
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
#' @returns Matrix of basis coefficients
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
  d <- ncol(data)
  if (is.null(d)){
    d <- 1
  }
  p <- ncol(basis)
  #ell <- nrow(basis)
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
    coeffs <- crossprod(y, x)
  }
  
  else if (attributes(weightinv)$diagonal == TRUE) {
    V <- t(basis) %*% (diag(weightinv) * basis)
    Q <- chol(V)
    y <- backsolve(Q, diag(p), transpose = TRUE)
    tmp <- t(basis) %*% (diag(weightinv) * data)
    x <- backsolve(Q, tmp, transpose = TRUE)
    coeffs <- crossprod(y, x)
  }
  
  else {
    V <- t(basis) %*% weightinv %*% basis
    Q <- chol(V)
    y <- backsolve(Q, diag(p), transpose = TRUE)
    x <- backsolve(Q, t(basis) %*% weightinv %*% data, transpose = TRUE)
    coeffs <- crossprod(y, x)
  }
  
  # Assign coefficient labels to columns
  coeffs <- t(coeffs)
  colnames(coeffs)[1:p] <- paste("C",1:p,sep="")
  
  return(coeffs)
}



#### see below - should add mean and variance? ####

#' Field reconstructions from coefficients
#'
#' Given a vector of coefficients for a basis, calculates the field
#'
#' @param coeffs Coefficient vector
#' @param basis Basis matrix
#' 
#' @returns Reconstructed field.
#'
#' @export
Recon <- function(coeffs, basis){
  
  if (is.null(ncol(basis))){
    q <- 1
  }
  else {
    q <- ncol(basis)
  }
  
  # Check whether provided with a single point or multiple
  if (is.null(ncol(coeffs))){
    # There's 2 possible cases: either q = 1 with n>= 1, or q > 1 and n = 1
    if (q == 1){
      n <- length(coeffs)
    }
    else {
      n <- 1
    }
  }
  else {
    n <- nrow(coeffs)
  }
  
  if (n > 1 & !(ncol(coeffs) == q)){
    stop('Inconsistency between number of coefficients and number of basis vectors')
  }
  
  else if (n == 1 & !(length(coeffs) == q)){
    stop('Inconsistency between number of coefficients and number of basis vectors')
  }
  
  if (q == 1 & n == 1){
    recons <- basis*as.numeric(coeffs)
  }
  else if (q > 1 & n == 1) {
    recons <- basis %*% as.numeric(coeffs)
  }
  else if (n > 1) {
    recons <- basis %*% t(coeffs)
  }
  
  return(recons)
}


# Recon <- function(coeff_expectation, coeff_variance, basis, ppe_mean){
#   q <- length(coeff_expectation)
#   n <- ncol(basis$v)
#   
#   field_exp <- ppe_mean + basis$v[,1:q] %*% coeff_expectation
#   field_var <- basis$v[,1:q] %*% diag(coeff_variance) %*% t(basis$v[,1:q])
#   
#   # Add variance due to removed (unemulated) basis vectors
#   deleted_var <- basis$d[-c(1:q,n)]^2 / (n-1) # gives the variance on each vector
#   deleted_bas <- basis$v[,-c(1:q,n)]
#   deleted_var <- deleted_bas %*% diag(deleted_var) %*% t(deleted_bas)
#   field_var <- field_var + deleted_var
#   
#   return(list(expectation = field_exp,
#               variance = field_var))
# }




# Change naming, do we need to centre automatically?

#' Project and reconstruct a given field
#'
#' Gives the reconstruction of a field using a basis, by projecting and back-projecting on this basis.
#'
#' @param field Vector over original field
#' @param basis Basis matrix
#'
#' @returns Reconstruction of the original field.
#'
#' @examples
#'
#' @export
ReconField <- function(field, basis, ...){
  coeffs <- Project(field, basis, ...)
  recons <- Recon(coeffs, basis)
  return(recons)
}





#' Matrix projection
#'
#' Projects a variance matrix onto a given basis
#'
#' @param mat A square matrix to be projected onto the basis
#' @param basis The basis to project with
#' @param weightinv The inverse of positive definite matrix W. If NULL, uses the standard projection, otherwise projects in the norm given by W.
#'
#' @returns The projection of the original matrix on the basis.
#'
#' @export
ProjectVar <- function(mat, basis, weightinv = NULL){
  if (is.null(weightinv)){
    proj <- t(basis) %*% mat %*% basis
  }
  else {
    V <- t(basis) %*% weightinv %*% basis
    Q <- chol(V)
    y <- backsolve(Q, diag(ncol(basis)), transpose = TRUE)
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
#' @returns The reconstruction error
#'
#' @export
ReconError <- function(obs, basis, weightinv = NULL, scale = FALSE){
  if (is.null(weightinv)){
    weightinv <- 0
    attr(weightinv, 'diagonal') <- attr(weightinv, 'identity') <- TRUE
  }
  field <- ReconField(obs, basis, weightinv)
  A <- c(obs) - field
  
  if (scale == TRUE){
    s <- length(c(obs))
  }
  else {
    s <- 1
  }
  if (attributes(weightinv)$diagonal == FALSE){
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


# Do we need a flag for centred?

#' Reconstruction errors
#' 
#' Calculates reconstruction error of observations for each truncation of a basis
#' 
#' @param obs (Centred) observation vector
#' @param basis The basis to project onto
#' 
#' @returns Vector of reconstruction errors for the observations
#' 
#' @export
errors <- function(basis, obs, weightinv=NULL){
  p <- ncol(basis)
  ell <- nrow(basis)
  err <- numeric(p)
  if (is.null(weightinv)){
    weightinv <- diag(ell)
    attr(weightinv, 'diagonal') <- attr(weightinv, 'identity') <- TRUE
  }
  for (i in 1:p){
    err[i] <- ReconError(obs, basis[,1:i], weightinv, scale = TRUE)
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
#' @returns The proportion of variability in the data that is explained by the basis
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
#' @returns The number of basis vectors required to explain vtot of the data.
#'
#' @export
ExplainT <- function(DataBasis, vtot = 0.95, weightinv = NULL){
  v <- 0
  q <- 0
  while (v < vtot & q < ncol(DataBasis$tBasis)){
    v <- VarExplained(DataBasis$tBasis[,1:(q+1)], DataBasis$CentredField, weightinv)
    q <- q + 1
  }
  return(q)
}






