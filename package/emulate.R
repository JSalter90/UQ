#' Functions for building and validating emulators

# Some preamble
n_cores <- ceiling(parallel::detectCores()/2)

cols_validate <- c('darkgrey', viridis::viridis(100)[31], viridis::viridis(100)[81])

cols_good <- viridis::viridis(100)[65]

#' Training emulators
#' 
#' A function for fitting 1 or more emulators.
#'
#' The default behaviour is to fit a single emulator, with output given by the final column of `tData`, and all preceding columns assumed as inputs.
#' 
#' The exception is if a column named C1 is found. If this exists, then this output and all subsequent columns are emulated, using all columns prior to C1 as inputs.
#' 
#' This default behaviour can be overriden by explicitly providing which columns relate to inputs and/or outputs.
#' 
#' @param tData A data frame containing the training data, with column(s) relating to all inputs and all outputs to be emulated. May also contain a noise column. May be the output of `ProcessData`.
#' @param method Which method to use for fitting the emulator(s). Defaults to rgasp `RobustGaSP`. Other options hetGP, gp, dgp (both use `dgpsi`).
#' @param input Column indices corresponding to emulator inputs. Defaults to `NULL`, in which case all columns prior to the first output column are assumed to be inputs.
#' @param output Column names corresponding to emulator outputs. Defaults to `NULL`, in which case either a column named C1 and all subsequent columns are emulated, or only the final column of `tData` is emulated.
#' @param output_cols Column indices corresponding to emulator outputs. Defaults to `NULL`, in which case behaves as `output`.
#' @param covariance Covariance function. Defaults to `matern5_2`.
#' @param nugget Logical, whether to estimate a nugget. Defaults to `TRUE`.
#' @param ... arguments for `BuildGasp`, `BuildHet`, `BuildGP` or `BuildDGP`.
#' 
#' @returns A list of emulators. See [BuildGasp()], [BuildHet()], [BuildGP()], [BuildDGP()] for the specific output given by each emulator type.
#' 
#' @export
FitEmulators <- function(tData, method = 'rgasp', input = NULL, output = NULL, output_cols = NULL, 
                         covariance = 'matern5_2', nugget = TRUE, ...){
  
  n <- nrow(tData)
  p1 <- ncol(tData)
  
  # If provided both output and output_cols, check these are consistent
  if (!(is.null(output)) & !(is.null(output_cols))){
    
    if (!(length(output) == length(output_cols))){
      stop('output and output_cols have different length, please provide consistent names/indices, or provide only one of these')
    }
    
    else {
      output1 <- sort(output)
      output2 <- sort(colnames(tData)[output_cols])
      if (!(all(output1 == output2))){
        stop('Inconsistent outputs selected by output and output_cols, please provide consistent names/indices, or provide only one of these')
      }
      rm(output1, output2)
    }
  }
  
  # If nothing specific provided, first search for column named C1, and emulate this and all subsequent columns
  if (is.null(output_cols) & is.null(output)){
    c1_ind <- which(colnames(tData) == 'C1')
    if (length(c1_ind) > 0){
      output_cols <- c1_ind:p1
      output <- colnames(tData)[output_cols]
    }
    
    # Otherwise select final column of tData to train emulator
    else {
      output_cols <- p1
      output <- colnames(tData)[output_cols]
    }
  }
  
  if (is.null(output_cols) & !is.null(output)){
    output_cols <- which(colnames(tData) %in% output)
  }
  
  if (!is.null(output_cols) & is.null(output)){
    output <- colnames(tData)[output_cols]
  }
  
  n_em <- length(output_cols) # number of emulators to be fitted
  
  # n_em may be empty if provided names are not contained in tData
  if (n_em < 1){
    stop('Please provide valid column indices or column names to be emulated')
  }
  
  # Attempt to select columns with too high index
  if (any(output_cols > p1)){
    stop('At least one provided index is higher than the number of dimensions of tData, please remove')
  }
  
  # If nothing specific provided, set inputs as all but emulated columns
  if (is.null(input)){
    input <- c(1:p1)[-output_cols]
  }
  
  if (length(input) == 0){
    stop('Please provide a data frame with more columns than those to be emulated')
  }
  
  # Filter down to required inputs and outputs
  tData <- tData[,c(input, output_cols)]
  
  # Relabel
  input <- 1:length(input)
  
  # Fit emulator(s)
  if (method %in% c('rgasp', 'gasp', 'Rgasp', 'Gasp')){
    if (n_em == 1){
      ems <- BuildGasp(Response = output[1], tData = tData, input = input, covariance = covariance, nugget = nugget, ...)
    }
    
    else {
      ems <- parallel::mclapply(output, 
                                function(k) BuildGasp(Response = k, tData = tData, input = input, covariance = covariance, nugget = nugget, ...),
                                mc.cores = n_cores)
    }
  }
  
  else if (method %in% c('gp')){
    if (n_em == 1){
      ems <- BuildGP(Response = output[1], tData = tData, input = input, covariance = covariance, nugget = nugget, ...)
    }
    
    else {
      ems <- lapply(output, 
                    function(k) BuildGP(Response = k, tData = tData, input = input, covariance = covariance, nugget = nugget, ...))
    }
  }
  
  else if (method %in% c('dgp', 'dgpsi')){
    if (n_em == 1){
      ems <- BuildDGP(Response = output[1], tData = tData, input = input, covariance = covariance, nugget = nugget, ...)
    }
    
    else {
      ems <- lapply(output, 
                    function(k) BuildDGP(Response = k, tData = tData, input = input, covariance = covariance, nugget = nugget, ...))
    }
  }
  
  else if (method %in% c('het', 'Het', 'hetgp', 'hetGP', 'HetGP')){
    if (n_em == 1){
      ems <- BuildHet(Response = output[1], tData = tData, input = input, covariance = covariance, ...)
    }
    
    else {
      ems <- parallel::mclapply(output, 
                                function(k) BuildHet(Response = k, tData = tData, input = input, covariance = covariance, ...),
                                mc.cores = n_cores)
    }
  }
  
  else {
    stop("method not recognised, please provide a valid selection (see ?FitEmulator)")
  }
  
  warning("Currently ignoring noise, fix this!!!")
  
  return(ems)
}


#' Replacing an emulator in a list
#' 
#' Retrains and replaces a specific emulator in a list of emulators (e.g., with a different emulator method or assumptions used).
#'
#' By default, uses the same set of inputs as the original emulator (as only these are stored in the training data). To override this, a different training dataset can be provided here.
#'
#' @param emulators A list of emulators, likely the output of [FitEmulators()].
#' @param Response The name of the output to be emulated.
#' @param tData A data frame containing the training data, with column(s) relating to all inputs and a column for the output to be emulated. Defaults to `NULL`, in which case inherits the training data from the previously fitted emulator. Can provide different set of inputs here.
#' @param method Which method to use for fitting the emulator(s). Defaults to rgasp `RobustGaSP`. Other options hetGP, gp, dgp (both use `dgpsi`).
#' @param ... Other options to pass to [FitEmulators()].
#' 
#' @returns A list of emulators, with the emulator relating to `Response` replaced with the new fitted emulator.
#' 
#' @export
RefitEmulator <- function(emulators, Response, tData = NULL, method = 'rgasp', ...){
  
  k <- length(emulators)
  all_responses <- unlist(lapply(1:k, function(k) emulators[[k]]$response))
  em_ind <- which(all_responses == Response)
  
  if (length(em_ind) == 0){
    stop('Response is not consistent with any component of emulators')
  }
  
  # Extract training data, if an alternative tData is not provided
  if (is.null(tData)){
    tData <- emulators[[em_ind]]$train_data
  }
  
  # Fit new emulator
  new_em <- FitEmulators(tData, method = method, ...)
  
  # Replace in emulators
  emulators[[em_ind]] <- new_em
  
  return(emulators)
}





#### DF vs matrix ####

#' Scale design
#' 
#' Scales input vectors to a common range.
#' 
#' @param design Design to be scaled, where each row as an input vector, each column is a particular input. The column names need to match with those provided in `InputRanges`.
#' @param InputRanges Matrix or data frame providing input ranges. Should contain 3 columns (name, minimum value, maximum value) which can have any names, with each row corresponding to a different input parameter.
#' @param range The range to scale inputs to. Defaults to (0,1).
#' 
#' @returns Scaled design, with same dimension as provided design.
#' 
ScaleInputs <- function(design, InputRanges, range = c(0,1)){
  
  # Rename columns
  InputRanges <- as.data.frame(InputRanges)
  colnames(InputRanges) <- c('Input', 'Min', 'Max')
  
  # Check all provided
  if (any(!(colnames(design) %in% InputRanges$Input))){
    stop('Please provide ranges for each input in design')
  }
  
  scaled_design <- design
  
  r1 <- range[1]
  r2 <- range[2]
  r <- r2 - r1
  
  for (i in 1:ncol(design)){
    ind <- which(InputRanges$Input == colnames(design)[i])
    scaled_design[,i] <- design[,i] - InputRanges$Min[ind]
    scaled_design[,i] <- scaled_design[,i] / ((InputRanges$Max[ind] - InputRanges$Min[ind])/r)
    scaled_design[,i] <- scaled_design[,i] + r1
  }
  
  return(scaled_design)
}


#' Unscale design
#' 
#' Unscales input vectors to their original range.
#' 
#' @param scaled_design The scaled design to be unscaled.
#' @param InputRanges Matrix or data frame providing input ranges. Should contain 3 columns (name, minimum value, maximum value) which can have any names, with each row corresponding to a different input parameter.
#' @param range The range the scaled inputs are on. Defaults to (0,1).
#' 
#' @returns Unscaled design, with same dimension as provided design.
#' 
UnscaleInputs <- function(scaled_design, InputRanges, range = c(0,1)){
  
  # Rename columns
  colnames(InputRanges) <- c('Input', 'Min', 'Max')
  
  # Check all provided
  if (any(!(colnames(scaled_design) %in% InputRanges$Input))){
    stop('Please provide ranges for each input in design')
  }
  
  design <- scaled_design
  
  r1 <- range[1]
  r2 <- range[2]
  r <- r2 - r1
  
  for (i in 1:ncol(scaled_design)){
    ind <- which(InputRanges$Input == colnames(scaled_design)[i])
    design[,i] <- scaled_design[,i] - r1
    design[,i] <- design[,i] * ((InputRanges$Max[ind] - InputRanges$Min[ind])/r)
    design[,i] <- design[,i] + InputRanges$Min[ind]
  }
  
  return(design)
}




#' Train/test sets
#' 
#' Splits a dataset into training and validation data randomly using a given training proportion.
#'
#' @param data Dataset to split into training and test sets.
#' @param train Proportion to assign to the training set. Defaults to 0.75.
#' @param seed Numeric seed that can be provided for reproducibility. Defaults to `NULL`.
#' 
#' @returns \item{train_data}{Training dataset.}
#' \item{train_inds}{Which rows of `data` were assigned to the training dataset.}
#' \item{val_data}{Test/validation dataset.}
#' \item{val_inds}{Which rows of `data` were assigned to the test/validation dataset.}
#' 
#' @export
TrainTestSplit <- function(data, train = 0.75, seed = NULL){
  n <- nrow(data)
  n_train <- ceiling(train * n)
  if (!(is.null(seed))){
    set.seed(seed)
  }
  samp_inds <- sample(1:n, size = n)
  train_inds <- samp_inds[1:n_train]
  val_inds <- samp_inds[-c(1:n_train)]
  
  train_data <- data[train_inds,]
  val_data <- data[val_inds,]
  
  return(list(train_data = train_data,
              train_inds = train_inds,
              val_data = val_data,
              val_inds = val_inds))
}







#' Building a single hetGP emulator
#'
#' Given tData object, fit an emulator for the selected coefficient
#'
#' @param Response a string indicating which output is being emulated. Must be consistent with a column of tData
#' @param tData data frame containing (in order): a) the design, b) a column containing noise, c) the basis coefficients
#' @param ... other arguments to pass to `mleHetGP`
#'
#' @returns \item{em}{An HetGP emulator}
#' \item{type}{Label indicating that the emulator was fitted with HetGP}
#' \item{train_data}{The subset of the data that was used to fit the emulator}
#' \item{validation_data}{The subset of the data that was not used. If training_prop = 1, this is empty}
#'
#' @export
BuildHet <- function(Response, tData, input = NULL, covariance = 'Matern5_2', ...){
  
  #### NOISE MIGHT NOT EXIST, NEED TO SELECT INPUTS BETTER ####
  # lastCand <- which(names(tData)=="Noise")
  # n <- dim(tData)[1]
  
  ind_response <- which(colnames(tData) == Response)
  
  # If inputs not provided explicitly, use all columns before Response
  if (is.null(input)){
    input <- 1:(ind_response-1)
  }
  
  tData_input <- tData[,input]
  tData_response <- tData[,ind_response]
  
  # Convert to correct naming for covariance
  if (covariance %in% c('matern5_2', 'matern_5_2', 'matern52', 'matern2.5', 'Matern5_2')){
    covariance <- 'Matern5_2'
  }
  
  else if (covariance %in% c('matern3_2', 'matern_3_2', 'matern32', 'matern1.5', 'Matern3_2')){
    covariance <- 'Matern3_2'
  }
  
  else if (covariance %in% c('sexp', 'squared_exp', 'pow_exp', 'gaussian', 'Gaussian')){
    covariance <- 'Gaussian'
  }
  
  else {
    stop('Please provide a valid covariance function for hetGP (Matern5_2, Matern3_2, Gaussian)')
  }
  
  # Process tData to handle replicates
  het_input <- hetGP::find_reps(X = as.matrix(tData_input), Z = tData_response)
  
  # Fit emulator
  em <- hetGP::mleHetGP(X = list(X0 = het_input$X0, Z0 = het_input$Z0, mult = het_input$mult),
                        Z = het_input$Z, covtype = covariance, ...)

  #train_data <- cbind(train_input, train_response)
  #colnames(train_data)[dim(train_data)[2]] <- Response
  #validation_data <- cbind(validation_input, validation_response)
  
  train_data <- tData[,c(input, ind_response)]
  
  return(list(em = em,
              method = 'het',
              response = Response,
              train_data = train_data))
}





#' Building a single GaSP emulator
#' 
#' Given tData object, fit an emulator for the selected coefficient
#' 
#' @param Response a string indicating which output is being emulated. Must be consistent with a column of tData
#' @param tData data frame containing (in order): a) the design, b) a column containing noise, c) the basis coefficients
#' @param mean_fn the structure allowed in the mean function. If NULL, fits a constant mean model. If 'linear', fits a model with a linear term in each of the inputs. If 'step', selects a mean function via stepwise regression
#' @param training_prop proportion of the data to use to fit the model, sampled at random
#' @param Fouriers if fitting a mean function via step, should Fourier terms be considered?
#' @param linModel defaults to NULL. If not, gives an object of type 'lm' to be used as the mean function
#' @param nugget_fit should a nugget be estimated? Defaults to TRUE
#' @param maxdf maximum number of terms allowed in the mean function, if step used. Defaults to 0.1*size of training data
#' 
#' @returns \item{em}{An rgasp emulator}
#' \item{em_lm}{If mean_fn = 'step', the regression model that was fitted}
#' \item{active}{If mean_fn = 'step', the variables that were deemed to be active}
#' \item{type}{Label indicating that the emulator was fitted with rgasp}
#' \item{mean_fn}{Label indicating the type of mean function used}
#' \item{train_data}{The subset of the data that was used to fit the emulator}
#' \item{validation_data}{The subset of the data that was not used. If training_prop = 1, this is empty}
#' 
#' @export
BuildGasp <- function(Response, tData, input = NULL, mean_fn = NULL, linModel = NULL, covariance = 'matern_5_2', nugget = TRUE, maxdf = NULL, ...){
  
  #### NOISE MIGHT NOT EXIST, NEED TO SELECT INPUTS BETTER ####
  #lastCand <- which(names(tData)=="Noise")
  n <- nrow(tData)
  ind_response <- which(colnames(tData) == Response)
  
  # If inputs not provided explicitly, use all columns before Response
  if (is.null(input)){
    input <- 1:(ind_response-1)
  }
  
  tData_input <- tData[,input]
  tData_response <- tData[,ind_response]
  train_data <- tData[,c(input, ind_response)]
  
  # Convert to correct naming for covariance
  if (covariance %in% c('matern5_2', 'matern_5_2', 'matern52', 'matern2.5', 'Matern5_2')){
    covariance <- 'matern_5_2'
  }
  
  else if (covariance %in% c('matern3_2', 'matern_3_2', 'matern32', 'matern1.5', 'Matern3_2')){
    covariance <- 'matern_3_2'
  }
  
  else if (covariance %in% c('sexp', 'squared_exp', 'pow_exp', 'gaussian', 'Gaussian')){
    covariance <- 'pow_exp'
  }
  
  else {
    stop('Please provide a valid covariance function for rgasp (matern_5_2, matern_3_2, pow_exp)')
  }
  
  # 'Trend' option defines mean function
  # e.g. have column of 1s for intercept, thereafter parameters
  # So for linear mean, just repeat design after column of 1s
  if (is.null(mean_fn)){
    em <- RobustGaSP::rgasp(design = tData_input, response = tData_response, kernel_type = covariance, nugget.est = nugget, ...)
  }
  
  else if (mean_fn == 'linear'){
    X <- cbind(rep(1,n), tData_input)
    em <- RobustGaSP::rgasp(design = tData_input, response = tData_response, trend = as.matrix(X), kernel_type = covariance, nugget.est = nugget, ...)
  }

  else if (mean_fn == 'lm'){
    em_lm <- mean_fn
    #X <- model.matrix(linModel)
    # In case linear model wasn't fitted with exact same data:
    tt <- stats::terms(linModel)
    tt <- stats::delete.response(tt)
    mm <- stats::model.frame(tt, tData_input)
    X <- stats::model.matrix(tt, mm)
    em <- RobustGaSP::rgasp(design = tData_input, response = tData_response, trend = as.matrix(X), kernel_type = covariance, nugget.est = nugget, ...)
  }
  
  if (!is.null(mean_fn)){
    if (mean_fn == 'lm'){
      return(list(em = em, lm = linModel, method = 'rgasp', response = Response, mean_fn = mean_fn, train_data = train_data))
    }
    else {
      return(list(em = em, method = 'rgasp', response = Response, mean_fn = mean_fn, train_data = train_data))
    }
  }
  if (is.null(mean_fn)){
    return(list(em = em, method = 'rgasp', response = Response, mean_fn = mean_fn, train_data = train_data))
  }
}






#' Building a single GP emulator using dgpsi
#' 
#' 
#'
#' @param Response 
#' @param tData 
#' @param input 
#' @param covariance 
#' @param nugget 
#' @param ... 
#'
#' @returns
#' @export
#'
#' @examples
BuildGP <- function(Response, tData, input = NULL, covariance = 'matern5_2', nugget = TRUE, ...){
  
  ind_response <- which(colnames(tData) == Response)
  
  # If input not provided explicitly, use all columns before Response
  if (is.null(input)){
    input <- 1:(ind_response-1)
  }
  
  tData_input <- tData[,input]
  tData_response <- tData[,ind_response]
  
  # Convert to correct naming for covariance
  if (covariance %in% c('matern5_2', 'matern_5_2', 'matern52', 'matern2.5', 'Matern5_2')){
    covariance <- 'matern2.5'
  }
  
  else if (covariance %in% c('sexp', 'squared_exp', 'pow_exp', 'gaussian', 'Gaussian')){
    covariance <- 'sexp'
  }
  
  else {
    stop('Please provide a valid covariance function for dgspi::gp (sexp, matern2.5)')
  }

  # Fit emulator
  em <- dgpsi::gp(as.matrix(tData_input), tData_response, name = covariance, nugget = 1e-8, nugget_est = nugget, verb = FALSE, ...)

  train_data <- tData[,c(input, ind_response)]
  
  return(list(em = em,
              method = 'gp',
              response = Response,
              train_data = train_data))
}


#' Building a single DGP emulator using dgpsi
#' 
#' 
#'
#' @param Response 
#' @param tData 
#' @param input 
#' @param covariance 
#' @param ... 
#'
#' @returns description
#' @export
#'
#' @examples
BuildDGP <- function(Response, tData, input = NULL, covariance = 'matern5_2', ...){
  
  ind_response <- which(colnames(tData) == Response)
  
  # If input not provided explicitly, use all columns before Response
  if (is.null(input)){
    input <- 1:(ind_response-1)
  }
  
  tData_input <- tData[,input]
  tData_response <- tData[,ind_response]
  
  # Convert to correct naming for covariance
  if (covariance %in% c('matern5_2', 'matern_5_2', 'matern52', 'matern2.5', 'Matern5_2')){
    covariance <- 'matern2.5'
  }
  
  else if (covariance %in% c('sexp', 'squared_exp', 'pow_exp', 'gaussian', 'Gaussian')){
    covariance <- 'sexp'
  }
  
  else {
    stop('Please provide a valid covariance function for dgspi::dgp (sexp, matern2.5)')
  }
  
  # Fit emulator
  em <- dgpsi::dgp(as.matrix(tData_input), tData_response, name = covariance, verb = FALSE, ...)
  
  train_data <- tData[,c(input, ind_response)]
  
  return(list(em = em,
              method = 'dgp',
              response = Response,
              train_data = train_data))
}



#' Wrapper of 1D predict function for multiple emulators
#' 
#' 
#' @param emulator 
#' @param design 
#' 
#' @returns description
#' 
#' @export
Predict <- function(emulator, design){

  # Find number of emulators
  if (!(is.null(emulator$method))){
    n_em <- 1
  }
  
  else {
    n_em <- length(emulator)
  }
  
  Expectation <- Variance <- matrix(0, nrow = nrow(design), ncol = n_em)
  
  if (n_em == 1){
    preds <- PredictSingle(emulator, design)
    Expectation[,1] <- preds$Expectation
    Variance[,1] <- preds$Variance
  }

  else {
    preds <- parallel::mclapply(1:n_em, 
                                function(k) PredictSingle(emulator[[k]], design),
                                mc.cores = n_cores)

    for (j in 1:n_em){
      Expectation[,j] <- preds[[j]]$Expectation
      Variance[,j] <- preds[[j]]$Variance
    }
  }
  
  return(list(Expectation = Expectation, 
              Variance = Variance))
}





#' Title
#'
#' @param emulator 
#' @param design 
#'
#' @returns
#'
#' @examples
PredictSingle <- function(emulator, design){
  
  # Check whether all required inputs are contained in design
  if (any(!(emulator$inputs %in% colnames(design)))){
    # extract which are missing
    missing_inputs <- ...
    stop(paste('Required input ', missing_inputs, 'not found in design'))
  }
  
  
  if (emulator$method == 'rgasp'){
    preds <- PredictGasp(emulator, design)
    preds$Expectation <- preds$mean
    preds$Variance <- preds$sd^2
  }
  
  else if (emulator$method %in% c('gp', 'dgp')){
    preds <- PredictDGP(emulator, design)
    preds$Expectation <- c(preds$mean)
    preds$Variance <- c(preds$var)
  }

  else if (emulator$method == 'het'){
    preds <- PredictHet(emulator, design)
    preds$Expectation <- preds$mean
    preds$Variance <- preds$sd2 + preds$nugs
  }
  
  return(preds)
}





#' Evaluating hetGP emulator predictions
#'
#' Given an object output by BuildHet, makes predictions for a set of inputs
#'
#' @param emulator an object output by BuildHet
#' @param design a data frame containing the input parameters, where each row is a point at which to evaluate the emulator
#'
#' @returns an object containing the mean, variance (with nugget removed), and nugget variance, at each input location
#'
#' @export
PredictHet <- function(emulator, design){
  # Ensure columns are ordered in the same way as when the emulator was trained
  col_names <- colnames(emulator$train_data)[-ncol(emulator$train_data)]
  design <- design[,col_names]
  
  try(hetGP::predict(emulator$em, as.matrix(design)), silent = TRUE)
  
  preds <- predict(emulator$em, as.matrix(design))
  return(preds)
}



#' Title
#'
#' @param emulator 
#' @param design 
#'
#' @returns
#' @export
#'
#' @examples
PredictDGP <- function(emulator, design){
  # Ensure columns are ordered in the same way as when the emulator was trained
  col_names <- colnames(emulator$train_data)[-ncol(emulator$train_data)]
  design <- design[,col_names]
  
  preds <- predict(emulator$em, as.matrix(design))$results
  return(preds)
}






#### need to handle noise ####

#' Evaluating GaSP emulator predictions
#' 
#' Given an object output by BuildGasp, makes predictions for a set of inputs
#' 
#' @param Design a data frame containing the input parameters, where each row is a point at which to evaluate the emulator
#' @param emulator an object output by BuildGasp (requires mean_fn, active etc.)
#' 
#' @returns an object containing the mean, standard deviation, and lower and upper bounds of the 95% posterior credible interval (see predict,rgasp-method)
#' 
#' @export
PredictGasp <- function(emulator, design){
  #### REMOVE, LIKELY REDUNDANT ####
  # if (class(emulator) == 'rgasp'){
  #   noise_ind <- colnames(design) == 'Noise'
  #   if (sum(noise_ind) > 0){
  #     noise_ind <- which(colnames(design) == 'Noise')
  #     design <- design[,c(1:(noise_ind-1))]
  #   }
  #   preds <- predict(emulator, design)
  # }
  
  if (is.null(emulator$mean_fn)){
    noise_ind <- colnames(design) == 'Noise'
    if (sum(noise_ind) > 0){
      noise_ind <- which(colnames(design) == 'Noise')
      design <- design[,c(1:(noise_ind-1))]
    }
    # Need to make sure columns are ordered in the same way as when the emulator was trained
    col_names <- colnames(emulator$train_data)[-ncol(emulator$train_data)]
    design <- design[,col_names]
    
    preds <- predict(emulator$em, design)
  }
  
  else if (emulator$mean_fn == 'linear'){
    noise_ind <- colnames(design) == 'Noise'
    if (sum(noise_ind) > 0){
      noise_ind <- which(colnames(design) == 'Noise')
      design <- design[,c(1:(noise_ind-1))]
    }
    # Need to make sure columns are ordered in the same way as when the emulator was trained
    col_names <- colnames(emulator$train_data)[-ncol(emulator$train_data)]
    design <- design[,col_names]
    X <- cbind(rep(1,dim(design)[1]), design)
    preds <- predict(emulator$em, design, testing_trend = as.matrix(X))
  }
  
  else if (emulator$mean_fn == 'lm'){
    tt <- stats::terms(emulator$lm)
    Terms <- stats::delete.response(tt)
    # Ensure columns are ordered in the same way as when the emulator was trained
    col_names <- colnames(emulator$train_data)[-ncol(emulator$train_data)]
    design <- design[,col_names]
    mm <- stats::model.frame(Terms, design, xlev = emulator$lm$xlevels)
    X <- stats::model.matrix(Terms, mm, contrasts.arg = emulator$lm$contrasts)
    preds <- predict(emulator$em, design, testing_trend = as.matrix(X))
  }
  return(preds)
}









#' Emulator validation
#' 
#' Predicting and plotting out-of-sample predictions vs the true output.
#'
#' @param emulator A single emulator, or a list of emulators.
#' @param valData Data frame of input vectors and true simulator output corresponding to each output in `emulator`. Set of input vectors used to evaluate the emulators, and plot the emulator predictions vs the truth.
#' @param interval Prediction interval to plot as error bars. Defaults to 0.95, can be in (0,1).
#' @param by_input 
#' @param by_index Whether the x axis should be the emulator prediction (`FALSE`, the default) or the training point index.
#' @param Obs If provided, adds horizontal and/or vertical dashed lines to show the location of the observation(s) relative to the emulator predictions. Defaults to `NULL`.
#'
#' @returns
#' @export
#'
#' @examples
Validate <- function(emulator, valData, interval = 0.95, by_input = FALSE, by_index = FALSE, Obs = NULL){
  
  if (interval <= 0 | interval >= 1){
    stop('interval must be between 0 and 1')
  }
  
  # Find number of emulators
  if (!(is.null(emulator$method))){
    n_em <- 1
  }
  
  else {
    n_em <- length(emulator)
  }
  
  if (n_em == 1){
    plot <- ValidateSingle(emulator, valData, interval, by_input, by_index, Obs = Obs)
  }
  
  else {
    plot <- parallel::mclapply(1:n_em, 
                               function(k) ValidateSingle(emulator[[k]], valData, interval, by_input, by_index, Obs = Obs),
                                mc.cores = n_cores)
  }
  
  if (n_em > 1 & !(by_input)){
    print(cowplot::plot_grid(plotlist = plot))
  }
  
  else if (n_em > 1 & by_input) {
    print(cowplot::plot_grid(plotlist = lapply(1:n_em, function(k) plot[[k]]$plot)))
  }
  
  return(plot)
}





#' Given a validation dataset, predicts and plots the mean and 95% uncertainty interval against the true output
#' 
#' @param emulator either a single output from BuildGasp, or a list of emulators
#' @param valData a validation data frame containing inputs and true output
#' @param IndivPars Create plots for each input
#'
#' @returns an object containing the mean, variance (with nugget removed), and nugget variance, at each input location
#'
ValidateSingle <- function(emulator, valData, interval = 0.95, by_input = FALSE, by_index = FALSE, Obs = NULL){
  
  if (interval <= 0 | interval >= 1){
    stop('interval must be between 0 and 1')
  }
  
  interval <- c((1-interval)/2, 1 - (1-interval)/2)
  
  valData <- valData[,colnames(emulator$train_data)]
  resp_ind <- ncol(valData)
  design <- valData[,-resp_ind]
  response <- valData[,resp_ind]
  
  # If only 1 input, make sure formatted properly
  if (ncol(emulator$train_data) == 2){
    design <- as.matrix(design, ncol = 1)
    colnames(design) <- colnames(emulator$train_data)[1]
  }
  
  if (emulator$method == 'rgasp'){
    preds <- PredictGasp(emulator, design)
    preds$lower95 <- preds$upper95 <- NULL
    preds$lower <- preds$mean + stats::qnorm(interval[1])*preds$sd
    preds$upper <- preds$mean + stats::qnorm(interval[2])*preds$sd
    preds$sd <- NULL
  }
  
  if (emulator$method %in% c('gp', 'dgp')){
    preds <- PredictDGP(emulator, design)
    preds$lower <- preds$mean + stats::qnorm(interval[1])*sqrt(preds$var)
    preds$upper <- preds$mean + stats::qnorm(interval[2])*sqrt(preds$var)
    preds$var <- preds$M <- NULL
  }
  
  if (emulator$method == 'het'){
    preds <- PredictHet(emulator, design)
    vars <- preds$sd2 + preds$nugs
    preds$lower <- preds$mean + stats::qnorm(interval[1])*sqrt(vars)
    preds$upper <- preds$mean + stats::qnorm(interval[2])*sqrt(vars)
    preds$sd2var <- preds$cov <- preds$sd2 <- preds$nugs <- NULL
  }
  
  upp <- max(c(preds$upper, response))
  low <- min(c(preds$lower, response))
  preds$truth <- response

  preds$In95 <- preds$truth >= preds$lower & preds$truth <= preds$upper
  perc_outside <- round(sum(preds$In95 == FALSE) / length(preds$In95) * 100, 1)
  
  cols <- cols_validate
  # Ensuring good points still coloured green if no points outside
  if (perc_outside == 0){
    cols[2:3] <- cols_validate[3]
  }
  
  if (by_index){
    preds$ind <- 1:length(preds$mean)
    plot1 <- ggplot(as.data.frame(preds), aes(x = .data$ind, y = .data$mean)) +
      geom_errorbar(aes(ymin = .data$lower, ymax = .data$upper), col = cols[1]) +
      geom_point(col = cols[1]) +
      geom_point(aes(x = .data$ind, y = .data$truth, col = .data$In95)) +
      scale_colour_manual(values = c(cols[2:3])) +
      labs(y = 'Prediction', x = 'Index', title = paste0('Outside ', 100*diff(interval), '% = ', perc_outside, '%')) +
      theme_bw() +
      theme(legend.position = 'none')
    
    if (!(is.null(Obs))){
      plot1 <- plot1 + geom_hline(aes(yintercept = as.numeric(Obs)), linetype = 'dashed')
    }
  }
  
  else {
    plot1 <- ggplot(as.data.frame(preds), aes(x = .data$truth, y = .data$mean, col = .data$In95)) +
      geom_errorbar(aes(ymin = .data$lower, ymax = .data$upper), col = cols[1]) +
      geom_point() +
      scale_colour_manual(values = c(cols[2:3])) +
      geom_abline(slope = 1, alpha = 0.6) +
      labs(y = 'Prediction', x = 'Truth', title = paste0('Outside ', 100*diff(interval), '% = ', perc_outside, '%')) +
      theme_bw() +
      theme(legend.position = 'none')
    
    if (!(is.null(Obs))){
      plot1 <- plot1 + 
        geom_hline(aes(yintercept = as.numeric(Obs)), linetype = 'dashed') + 
        geom_vline(aes(xintercept = as.numeric(Obs)), linetype = 'dashed')
    }
  }

  if (by_input == TRUE){
    plot_data <- data.frame(design, as.data.frame(preds))
    plot_data$ind <- NULL
    plot_data <- reshape2::melt(plot_data, id.vars = c('mean', 'lower', 'upper', 'truth', 'In95'))
    plot2 <- ggplot(plot_data, aes(x = .data$value, y = .data$mean)) +
      geom_errorbar(aes(ymin = .data$lower, ymax = .data$upper), col = cols[1]) +
      geom_point(col = cols[1]) +
      geom_point(aes(x = .data$value, y = .data$truth, col = .data$In95)) +
      facet_wrap(vars(.data$variable), scales = 'free_x') +
      scale_colour_manual(values = c(cols[2:3])) +
      labs(y = 'Prediction', x = 'Input', title = paste0('Outside ', 100*diff(interval), '% = ', perc_outside, '%')) +
      theme_bw() +
      theme(legend.position = 'none')
    
    if (!(is.null(Obs))){
      plot2 <- plot2 + geom_hline(aes(yintercept = as.numeric(Obs)), linetype = 'dashed')
    }
  }
  
  if (by_input == FALSE){
    return(plot1)
  }
  
  else{
    return(list(plot = plot1,
                input = plot2))
  }
}







#' Emulator leave-one-out predictions
#' 
#' Plotting emulator leave-one-out predictions vs the true output.
#'
#' @param emulator A single emulator, or a list of emulators.
#' @param interval Prediction interval to plot as error bars. Defaults to 0.95, can be in (0,1).
#' @param Obs If provided, adds horizontal and/or vertical dashed lines to show the location of the observation(s) relative to the emulator predictions. Defaults to `NULL`.
#' @param by_index Whether the x axis should be the emulator prediction (`FALSE`, the default) or the training point index.
#'
#' @returns Leave-one-out plots for each emulator.
#' @export
#'
#' @examples
LeaveOneOut <- function(emulator, interval = 0.95, Obs = NULL, by_index = FALSE){
  
  # Find number of emulators
  if (!(is.null(emulator$method))){
    n_em <- 1
  }
  
  else {
    n_em <- length(emulator)
  }
  
  if (n_em == 1){
    plot <- LeaveOneOutSingle(emulator, Obs = Obs, by_index = by_index)
  }
  
  else {
    plot <- parallel::mclapply(1:n_em, 
                               function(k) LeaveOneOutSingle(emulator[[k]], Obs = Obs[k], interval = interval, by_index = by_index),
                               mc.cores = n_cores)
  }
  
  if (n_em > 1){
    print(cowplot::plot_grid(plotlist = plot))
  }

  return(plot)
}






#' Leave-one-out
#' 
#' @param emulator 
#' 
#' @returns description
#' 
LeaveOneOutSingle <- function(emulator, interval = 0.95, by_index = FALSE, Obs = NULL){

  if (interval <= 0 | interval >= 1){
    stop('interval must be between 0 and 1')
  }
  
  interval <- c((1-interval)/2, 1 - (1-interval)/2)
  
  response <- emulator$train_data[,dim(emulator$train_data)[2]]
  
  if (emulator$method == 'rgasp'){
    loo_preds <- RobustGaSP::leave_one_out_rgasp(emulator$em)
    loo_preds$lower <- loo_preds$mean + stats::qnorm(interval[1])*loo_preds$sd
    loo_preds$upper <- loo_preds$mean + stats::qnorm(interval[2])*loo_preds$sd
  }
  
  if (emulator$method %in% c('gp', 'dgp')){
    loo_preds <- dgpsi::validate(emulator$em, verb = FALSE)$loo
    loo_preds$lower <- loo_preds$lower
    loo_preds$upper <- loo_preds$upper
  }
  
  if (emulator$method == 'het'){
    loo_preds <- hetGP::LOO_preds(emulator$em)
    loo_preds$lower <- loo_preds$mean + stats::qnorm(interval[1])*sqrt(loo_preds$sd2)
    loo_preds$upper <- loo_preds$mean + stats::qnorm(interval[2])*sqrt(loo_preds$sd2)
  }

  loo_preds$truth <- response
  upp <- max(c(loo_preds$upper, response))
  low <- min(c(loo_preds$lower, response))
  
  loo_preds$In95 <- loo_preds$truth >= loo_preds$lower & loo_preds$truth <= loo_preds$upper
  perc_outside <- round(sum(loo_preds$In95 == FALSE) / length(loo_preds$In95) * 100, 1)
  
  cols <- cols_validate
  # Ensuring good points still coloured green if no points outside
  if (perc_outside == 0){
    cols[2:3] <- cols_validate[3]
  }
  
  if (by_index){
    loo_preds$ind <- 1:length(loo_preds$mean)
    plot <- ggplot(as.data.frame(loo_preds), aes(x = .data$ind, y = .data$mean)) +
      geom_errorbar(aes(ymin = .data$lower, ymax = .data$upper), col = cols[1]) +
      geom_point(col = cols[1]) +
      geom_point(aes(x = .data$ind, y = .data$truth, col = .data$In95)) +
      scale_colour_manual(values = c(cols[2:3])) +
      labs(y = 'Prediction', x = 'Input', title = paste0('Outside ', 100*diff(interval), '% = ', perc_outside, '%')) +
      theme_bw() +
      theme(legend.position = 'none')
    
    if (!(is.null(Obs))){
      plot <- plot + geom_hline(aes(yintercept = as.numeric(Obs)), linetype = 'dashed')
    }
  }
  
  else {
    plot <- ggplot(as.data.frame(loo_preds), aes(x = .data$truth, y = .data$mean, col = .data$In95)) +
      geom_errorbar(aes(ymin = .data$lower, ymax = .data$upper), col = cols[1]) +
      geom_point() +
      scale_colour_manual(values = c(cols[2:3])) +
      geom_abline(slope = 1, alpha = 0.6) +
      labs(y = 'Prediction', x = 'Truth', title = paste0('Outside ', 100*diff(interval), '% = ', perc_outside, '%')) +
      theme_bw() +
      theme(legend.position = 'none')
    
    if (!(is.null(Obs))){
      plot <- plot + 
        geom_hline(aes(yintercept = as.numeric(Obs)), linetype = 'dashed') + 
        geom_vline(aes(xintercept = as.numeric(Obs)), linetype = 'dashed')
    }
  }

  return(plot)
}




#' Validating predictions based on aggregations
#'
#' For a set of samples, and chosen output indices/locations, compares the emulator expectation and variance for this subset to the true values
#'
#' @param Samples A set of samples from the basis emulator, reconstructed to the original field. Usually the output of `BasisEmSamples()`
#' @param Truth The corresponding true model outputs, for comparison. Must have an ordering of rows and columns consistent with `Samples` (in terms of input and output locations)
#' @param output_inds Indices of locations (across the output field) to average over. Defaults to NULL, which uses all outputs
#'
#' @returns Validation plot comparing truth and emulator samples
#'
#' @export
ValidateSummary <- function(emulator, design, DataBasis, interval = 0.95, output_inds = NULL, data_inds = NULL, plot_sum = FALSE){
  
  if (interval <= 0 | interval >= 1){
    stop('interval must be between 0 and 1')
  }
  
  interval <- c((1-interval)/2, 1 - (1-interval)/2)
  
  if (is.null(data_inds)){
    data_inds <- 1:nrow(design)
  }
  
  # Predict
  preds <- Predict(emulator, design[data_inds,])

  # Draw samples
  Samples <- BasisEmSamples(preds, DataBasis)
  
  # Construct truth
  Truth <- DataBasis$CentredField[,data_inds] + DataBasis$EnsembleMean
  
  ell <- dim(Samples)[1]
  ns <- dim(Samples)[2]
  
  # If not provided with a subset of locations, use all
  if (is.null(output_inds)){
    output_inds <- 1:ell
  }
  
  # Find average or sum across output_inds for each sample, run
  if (plot_sum){
    SamplesSummary <- apply(Samples[output_inds,,], c(2,3), sum)
  }
  
  else {
    SamplesSummary <- apply(Samples[output_inds,,], c(2,3), mean)
  }
  
  
  # Find mean/95% interval of averages for plotting
  SamplesMean <- apply(SamplesSummary, 2, mean)
  SamplesLower <- apply(SamplesSummary, 2, stats::quantile, prob = interval[1])
  SamplesUpper <- apply(SamplesSummary, 2, stats::quantile, prob = interval[2])
  
  if (!(nrow(Truth) == length(output_inds))){
    Truth <- Truth[output_inds,]
  }
  
  if (length(output_inds) == 1){
    TruthMean <- Truth
  }
  
  else {
    if (plot_sum){
      TruthMean <- apply(Truth, 2, sum)
    }
    else {
      TruthMean <- apply(Truth, 2, mean)
    }
  }
  
  plot_data <- data.frame(Mean = SamplesMean,
                          Lower = SamplesLower,
                          Upper = SamplesUpper,
                          Truth = TruthMean)
  
  plot_data$In95 <- plot_data$Truth >= plot_data$Lower & plot_data$Truth <= plot_data$Upper
  perc_outside <- round(sum(plot_data$In95 == FALSE) / length(plot_data$In95) * 100, 1)
  
  cols <- cols_validate
  # Ensuring good points still coloured green if no points outside
  if (perc_outside == 0){
    cols[2:3] <- cols_validate[3]
  }
  
  plot <- ggplot(plot_data, aes(.data$Truth, .data$Mean, col = .data$In95)) +
    geom_errorbar(aes(ymin = .data$Lower, ymax = .data$Upper), col = cols[1]) +
    geom_point() +
    scale_colour_manual(values = c(cols[2:3])) +
    ggplot2::geom_abline(slope = 1, alpha = 0.6) +
    labs(y = 'Prediction', title = paste0('Outside ', 100*diff(interval), '% = ', perc_outside, '%')) +
    theme_bw() +
    theme(legend.position = 'none')
  
  return(plot)
}




