#' Calculate a set of vectors from weights and a basis
#'
#' Give basis vectors from linear combinations
#'
#' @param weights A vector of weights
#' @param basis Original basis
#'
#' @return New basis vector(s)
#'
#' @export
ReconBasis <- function(weights, basis){
  n <- dim(basis)[2]
  q <- length(weights) / n
  if (q == 1){
    new.basis <- as.vector(tensor(basis, weights, 2, 1))
  }
  else {
    dim(weights) <- c(n, q)
    new.basis <- tensor(basis, weights, 2, 1)
  }
  return(new.basis)
}



#### More flexibility in specification, e.g. v, time allowed, make prior clearer, remove months etc. ####
#' Finding a calibration-optimal basis rotation
#'
#' Given a DataBasis object, observations, matrix W, vector v, applies the optimal rotation algorithm to find a basis more suitable for calibration.
#'
#' @param DataBasis An object containing lxn ensemble data and the lxn basis that will be rotated
#' @param obs A vector of length l with observations on the same scale as the ensemble
#' @param kmax Maximum number of iterations allowed (defaults to 5)
#' @param weightinv Inverse of positive definite weight matrix W. If set = NULL, uses the identity.
#' @param v Vector of minimum proportion of the ensemble data to be explained by the corresponding rotated basis vector
#' @param vtot Minimum proportion of ensemble variability to be explained by the truncated basis
#' @param MaxTime Maximum time (in seconds) to run the optimiser at each iteration.
#'
#' @return \item{tBasis}{Full rotated basis}
#' \item{CentredField}{The ensemble that was passed into the function}
#' \item{EnsembleMean}{The ensemble mean}
#' \item{scaling}{Initial scaling applied to the data}
#' \item{RW}{The reconstruction error after each iteration}
#' \item{VarExp}{The variance explained by each rotated basis vector}
#' \item{Opt}{Linear combination of the basis that gave rotated basis vectors} #### only true for first really ####
#'
#'@examples # First run an ensemble of idealised function fn
#'
#' n <- 60
#' sample <- as.data.frame(2*maximinLHS(n,6) - 1)
#' colnames(sample) <- c("x1","x2","x3","x4","x5","x6")
#' output <- array(c(rep(0,100*n)), dim=c(10,10,n))
#' for (i in 1:n){
#'   output[,,i] <- fn(as.numeric(sample[i,]))
#' }
#' dim(output) <- c(100, n)
#'
#' DataBasis <- MakeDataBasis(data = output, RemoveMean = TRUE)
#'
#' # Define the observations as a known value of x, plus some noise
#'
#' obs <- c(fn(c(0.7,0.01,0.01,0.25,0.8,-0.9)) + rnorm(100, mean = 0, sd = 0.1))
#' obs <- obs - DataBasis$EnsembleMean # centre observations so that comparable to the data
#'
#' # Look at the VarMSEplot for the SVD basis
#'
#' vSVD <- VarMSEplot(DataBasis = DataBasis, obs = obs)
#'
#' # Perform a rotation
#'
#' RotatedBasis <- RotateBasis(DataBasis = DataBasis, obs = obs, kmax = 3)
#'
#' # Editing the variance constraints so that the first vector explains at least 40% of the ensemble, at least 10% for later vectors
#'
#' RotatedBasis <- RotateBasis(DataBasis = DataBasis, obs = obs, kmax = 3, v = c(0.4,0.1,0.1))
#'
#' # So far assumed that W is the identity. Now add structure to W
#'
#'
#'
#' @export
RotateBasis <- function(DataBasis, obs, kmax = 5, weightinv = NULL, v = c(rep(0.1,5)), vtot = 0.95, prior = NULL,
                        StoppingRule = TRUE, MaxTime = 60, ...){
  data <- DataBasis$CentredField
  basis <- DataBasis$tBasis
  if (!dim(data)[1] == dim(basis)[1]){
    stop("Dimension of ensemble and basis (l) are not the same")
  }
  obs <- c(obs)
  if (!length(obs) == dim(basis)[1]){
    stop("Observations not the correct dimension (l)")
  }
  l <- dim(basis)[1]
  n <- dim(data)[2]
  minRw <- ReconError(obs, basis, weightinv)
  if (is.null(prior)){
    prior <- c(1:dim(basis)[2])
  }
  basis <- basis[,prior]
  mse <- var <- numeric(kmax)
  x <- NULL
  new.basis <- NULL
  if (is.null(weightinv)){
    var_sum <- crossprod(c(data))
  }
  else {
    if (attributes(weightinv)$diagonal == TRUE){
      var_sum <- sum(t(data)^2 %*% diag(weightinv))
    }
    else {
      var_sum <- sum(diag(t(data) %*% weightinv %*% data))
    }
  }
  if (is.null(DataBasis$Q)){
    Q <- NULL
    Lambda <- NULL
  }
  else {
    Q <- DataBasis$Q
    Lambda <- DataBasis$Lambda
  }
  for (i in 1:kmax){
    p <- dim(basis)[2]
    if (is.null(weightinv)){
      psi <- t(basis) %*% diag(l) %*% basis
    }
    else {
      psi <- t(basis) %*% weightinv %*% basis
    }
    opt <- GenSA(c(1, rep(0, p-1)), WeightOptim,  lower = rep(-1, p*1),
                 upper = rep(1, p*1), basis = basis, obs = obs, data = data, weightinv = weightinv,
                 v = v[i], newvectors = new.basis, total_sum = var_sum, psi = psi, control = list(max.time = MaxTime), ...)
    best.patterns <- cbind(new.basis, ReconBasis(opt$par, basis))
    rank <- min(n, l)
    basis <- ResidBasis(best.patterns, data, weightinv, Q, Lambda)[,1:rank]
    x <- c(x, opt$par)
    q <- ExplainT(DataBasis, vtot, weightinv)
    mse[i] <- ReconError(obs, basis[,1:q], weightinv)
    var[i] <- VarExplained(basis[,i], data, weightinv)
    new.basis <- cbind(new.basis, basis[,i])
    basis <- basis[,-(1:i)]
    if (round(mse[i],4) == round(minRw,4)) break #### INSERT STOPPING RULE
  }
  new.basis <- cbind(new.basis, basis)[,1:rank]
  return(list(tBasis = new.basis, CentredField = DataBasis$CentredField,
              EnsembleMean = DataBasis$EnsembleMean, scaling = DataBasis$scaling,
              RW = mse, VarExp = var, Opt = x))
}


#' Find the residual basis
#'
#' Given basis vectors and data, calculate the residual basis to complete the basis for the data
#'
#' @param basisvectors Basis vector(s) to project the data onto
#' @param data Data matrix for to find a basis for
#'
#' @return The full basis for the data, i.e. the basis vectors passed into the function, with the residual basis vectors appended
#'
#' @export
ResidBasis <- function(basisvectors, data, weightinv = NULL, ...){
  if (is.null(weightinv)){
    basisvectors <- orthonormalization(basisvectors,basis=FALSE,norm=TRUE)
  }
  else {
    newvector <- basisvectors[,dim(basisvectors)[2]]
    basisvectors[,dim(basisvectors)[2]] <- newvector / as.numeric(sqrt(t(newvector)%*%weightinv %*% newvector))
  }
  l <- dim(data)[1]
  n <- dim(data)[2]
  recons <- matrix(numeric(l*n), nrow=l)
  for (i in 1:n){
    recons[,i] <- ReconObs(data[,i],basisvectors, weightinv)
  }
  resids <- data - recons
  if (is.null(weightinv)){
    svd.resid <- svd(t(resids))
  }
  else {
    svd.resid <- wsvd(t(resids), weightinv = weightinv, ...)
  }
  new.basis <- cbind(basisvectors, svd.resid$v)[,1:min(l,n)]
  return(new.basis)
}



#' Produce a VarMSEplot
#'
#' Calculates and plots the reconstruction error and proportion of variability explained for each truncated basis
#'
#'  @param DataBasis A DataBasis object, containing the data and basis to be used to produce the VarMSEplot
#'  @param obs Observation vector
#'  @param RecVarData if given, the output given by a previous call of VarMSEplot, so that plots can be reproduced without rerunning potentially slow calculations
#'  @param weightinv Inverse of W used for calculating the reconstruction error. If NULL, calculates the mean squared error
#'  @param min.line If TRUE, plots a solid horizontal line at the minimum value of the reconstruction error
#'  @param bound If TRUE, plots a dotted horizontal line at the value of the history matching bound
#'  @param qmax if not NULL, value for the maximum size of the truncated basis to consider. Useful if W is non-diagonal, and either n or l is large
#'
#'  @return A VarMSEplot, and a matrix containing the plotted data, arranged in columns as (Reconstruction error, proportion of variability explained)
#'
#'  @export
VarMSEplot <- function(DataBasis, obs, RecVarData = NULL, weightinv=NULL, min.line=TRUE, bound=TRUE, qmax = NULL, ...){
  if (!is.null(RecVarData)){
    PlotData <- RecVarData
    qmax <- p <- dim(PlotData)[1]
  }
  else {
    p <- dim(DataBasis$tBasis)[2]
    if (is.null(qmax)){
      qmax <- p
    }
    PlotData <- matrix(numeric(qmax*2), nrow=qmax)
    if (!is.null(DataBasis$scaling)){
      PlotData[,1] <- errors(DataBasis$tBasis[,1:qmax], obs, weightinv)*DataBasis$scaling^2
    }
    else {
      PlotData[,1] <- errors(DataBasis$tBasis[,1:qmax], obs, weightinv)
    }
    if (is.null(weightinv)){
      var_sum <- crossprod(c(DataBasis$CentredField))
    }
    else {
      var_sum <- sum(diag(t(DataBasis$CentredField) %*% weightinv %*% DataBasis$CentredField))
    }
    for (i in 1:qmax){
      PlotData[i,2] <- VarExplained(DataBasis$tBasis[,1:i], DataBasis$CentredField, weightinv, total_sum = var_sum)
    }
  }
  if (qmax < p){
    FinalR <- ReconError(obs, DataBasis$tBasis, weightinv)
    PlotData <- rbind(PlotData, c(FinalR, 1))
    plotseq <- c(1:qmax, p)
  }
  else {
    plotseq <- 1:p
  }
  plot(plotseq, PlotData[,1], type="l", col="red",xlab = expression(k), ylab = '', ...)
  mtext(side = 2, line = 2.5, expression(paste("R "[bold(W)], " (", bold(B)[k], ",", bold(z), ")", " / l")), las = 3, cex = 0.8)
  if (min.line)
    abline(h = min(PlotData[,1]), col=alpha("black", 0.7), lty=4)
  if (bound == TRUE){
    abline(h = qchisq(0.995, length(obs))/length(obs), lty=2)
  }
  par(new = TRUE)
  plot(plotseq, PlotData[,2], type="l", axes=FALSE, xlab=NA, ylab=NA, col="blue", ylim=c(0,1))
  axis(side = 4)
  mtext(side = 4, line = 2.5, expression(paste("V(", bold(B)[k], ",", bold(F), ")")), las=3)
  return(PlotData)
}




#' Function used within optimiser for minimisation of the reconstruction error for rotated basis, subject to constraints
#'
#' Given a vector of weights, gives the new basis vector as a linear combination of the original basis, and calculates the reconstruction error, subject to a variability constraint
#'
#' @param x Vector giving the linear combination of the basis to use
#' @param basis The basis that is being rotated
#' @param obs Observation vector
#' @param data Ensemble data
#' @param weightinv Inverse of matrix W
#' @param v The proportion of variability to be explained by the basis vector
#' @param total_sum Common denominator used to calculate VarExplained
#' @param psi As new basis is linear combinarion of original, if pass psi = t(basis) %*% weightinv %*% basis adds efficiency
#' @param newvectors If the reconstruction error should account for any previous basis vectors
#'
#' @return The reconstruction error
#'
#' @export
WeightOptim <- function(x, basis, obs, data, weightinv, v = 0.1, total_sum = NULL, psi = NULL, newvectors = NULL){
  new.basis <- as.vector(tensor(basis, x, 2, 1))
  if (is.null(newvectors) == FALSE){
    new.basis <- cbind(newvectors, new.basis)
  }
  if (is.null(newvectors) == TRUE){
    v_new <- VarExplained(new.basis, data, weightinv, total_sum, psi, basis_lincom = x)
  }
  else {
    v_new <- VarExplained(new.basis[,dim(new.basis)[2]], data, weightinv, total_sum, psi, basis_lincom = x)
  }
  if (v_new < v){
    y <- 999999999
  }
  else {
    y <- ReconError(obs, new.basis, weightinv)
  }
  return(y)
}

#' Idealised function with spatial output
#'
#' A 6 parameter toy function that gives output over a 10x10 field.
#'
#' @param x A vector of 6 values, corresponding to x_1, ..., x_6.
#'
#' @return A 10x10 field.
#'
#' @examples n <- 60
#' sample <- as.data.frame(2*maximinLHS(n,6) - 1)
#' colnames(sample) <- c("x1","x2","x3","x4","x5","x6")
#' output <- array(c(rep(0,100*n)), dim=c(10,10,n))
#' for (i in 1:n){
#'   output[,,i] <- fn(as.numeric(sample[i,]))
#' }
#' dim(output) <- c(100, n) # vectorising spatial output, so each column of the data matrix is a single realisation of the function
#'
#' @export
fn <- function(x){
  basis1 <- basis2 <- basis3 <- basis4 <- basis5 <- basis6 <- basis7 <- basis8 <- matrix(c(rep(0,100)),nrow=10)
  for (i in 1:10){basis1[i,i] <- 1}
  basis1[10,1] <- basis1[10,2] <- basis1[9,1] <- basis1[1,9] <- basis1[1,10] <- basis1[2,10] <- -1
  basis1[8,1] <- basis1[9,2] <- basis1[10,3] <- basis1[1,8] <- basis1[2,9] <- basis1[3,10] <- -1
  for (i in 1:9){basis2[i,i+1] <- 1}
  for (i in 1:8){basis3[1,i] <- 1}
  for (i in 1:5){basis3[2,i] <- 1}
  for (i in 1:3){basis3[3,i] <- 1}
  basis3[9,9] <- basis3[10,10] <- basis3[9,10] <- basis3[10,9] <- basis3[8,9] <- basis3[8,10] <- basis3[7,2] <- basis3[7,3] <- basis3[7,4] <- basis3[8,2] <- basis3[8,3] <- basis3[8,4] <- -1
  basis3[4,1] <- basis3[5,2] <- basis3[5,3] <- basis3[5,4] <- basis3[6,2] <- basis3[6,3] <- basis3[6,4] <- 1
  for (i in 3:7){basis4[10,i] <- basis4[9,i] <- 1}
  basis4[10,3] <- 0
  basis4[6,1] <- basis4[6,2] <- basis4[7,1] <- basis4[6,3] <- basis4[1,1] <- basis4[1,2] <- basis4[2,1] <- -1
  basis4[5,2] <- basis4[5,3] <- basis4[2,3] <- basis4[2,4] <- basis4[3,3] <- 1
  for (i in 8:10){basis5[4,i] <- basis5[5,i] <- basis5[6,i] <- basis5[7,i] <- 1}
  for (i in 6:8){basis5[3,i] <- basis5[2,i] <- basis5[1,i] <- -1}
  basis5[4,7] <- basis5[4,6] <- -1
  basis5[1,8] <- 0
  for (i in 8:9){basis5[i,10] <- -1}
  for (i in 6:8){basis6[i,5] <- basis6[i,6] <- 1}
  basis6[5,6] <- basis6[5,7] <- basis6[8,5] <- basis6[8,4] <- 1
  for (i in 6:7){basis6[i,2] <- basis6[i,3] <- -1}
  for (i in 6:8){basis6[10,i] <- -1}
  basis6[9,9] <- basis6[9,10] <- basis6[1,4] <- basis6[1,5] <- basis6[8,2] <- basis6[6,9] <- -1
  basis7[2,10] <- basis7[3,10] <- basis7[3,9] <- basis7[3,8] <- basis7[10,8] <- basis7[9,8] <- basis7[8,8] <- basis7[9,9] <- basis7[7,6] <- basis7[7,9] <- basis7[5,7] <- basis7[3,2] <- basis7[4,2] <- basis7[1,8] <- basis7[2,8] <- basis7[2,2] <- basis7[8,3] <- basis7[7,3] <- basis7[10,4] <- basis7[9,4] <- basis7[8,7] <- 1
  basis7[7,5] <- basis7[6,5] <- basis7[6,6] <- basis7[7,4] <- basis7[8,1] <- basis7[8,2] <- basis7[1,3] <- basis7[1,4] <- basis7[10,6] <- basis7[10,5] <- basis7[6,9] <- basis7[6,10] <- basis7[1,6] <- basis7[2,6] <- basis7[2,7] <- basis7[10,9] <- -1
  basis8[3,5] <- basis8[3,6] <- basis8[4,5] <- basis8[4,4] <- basis8[5,10] <- basis8[6,10] <- basis8[1,9] <- basis8[9,7] <- 1
  basis8[9,5] <- basis8[8,5] <- basis8[7,5] <- basis8[10,7] <- basis8[7,2] <- basis8[5,4] <- basis8[7,1] <- basis8[10,1] <- basis8[7,10] <- basis8[7,7] <- basis8[5,1] <- basis8[4,3] <- basis8[6,7] <- -1
  fx <- 3*(10*x[2]^2*basis2 + 5*x[3]^2*basis2 + (x[3] + 1.5*x[1]*x[2])*basis3 + 2*x[2]*basis4 + x[3]*x[1]*basis5 +
             (x[2]*x[1])*basis6 + (x[2]^3)*basis7 + ((x[2] + x[3])^2)*basis8 + 2) + 1.5*dnorm(x[4], 0.2, 0.1)*basis1*(x[5]/(1.3+x[6]))
  noise <- matrix(rnorm(100,0,0.05),nrow=10)
  fx <- fx + noise
  return(fx)
}





