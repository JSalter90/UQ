# Some additional basis functionality, that will be combined with other source code at some point
# Mostly plotting, but also sampling from emulator posteriors
library(ggplot2)
library(reshape2)

# Colours used in several validation functions
my.cols <- c('darkgrey', viridis::viridis(100)[c(31,81)])

#' Plotting emulator samples
#'
#' Takes samples (or a summary of these) from BasisEmSamples and plots, adds the truth if this is provided
#'
#' @param Samples 
#' @param inds 
#' @param Truth 
#' @param input_values 
#' @param input_name 
#' @param output_name 
#' 
#' @return 
#' @import ggplot2
#' 
#' @export
PlotSamples <- function(Samples, inds = NULL, Truth = NULL, input_values = NULL, input_name = NULL, output_name = NULL){
  
  if (is.null(input_name)){
    input_name <- 'Input'
  }
  
  if (is.null(output_name)){
    output_name <- 'Output'
  }

  # If all samples are provided
  if (!(is.null(dim(Samples)))){
    
    ell <- dim(Samples)[1]
    ns <- dim(Samples)[2]
    runs <- dim(Samples)[3]
    
    # If only have samples from a single input
    if (is.na(runs)){
      inds <- 1
      runs <- 1
    }
    
    # If don't provide which inputs to plot, and there's more than 1 option, set as all
    if (is.null(inds)){
      inds <- 1:runs
    }
    
    # If have multiple inputs and a set of indices was provided
    if (runs > 1 & !(is.null(inds))){
      Samples <- Samples[,,inds]
    }

    if (is.null(input_values)){
      input_values <- 1:ell
    }

    plot_data <- data.frame(Input = input_values, 
                            Output = c(Samples), 
                            s = rep(1:ns, each = ell),
                            Run = rep(inds, each = ell*ns))
    
    plot <- ggplot(plot_data, aes(.data$Input, .data$Output, col = as.factor(.data$s))) + 
      geom_line() +
      geom_line(data = data.frame(aggregate(Output ~ Input + Run, plot_data, mean)), col = 'red', size = 1.25) +
      geom_line(data = data.frame(aggregate(Output ~ Input + Run, plot_data, quantile, probs = 0.025)), col = 'red', size = 1.25, linetype = 'dashed') +
      geom_line(data = data.frame(aggregate(Output ~ Input + Run, plot_data, quantile, probs = 0.975)), col = 'red', size = 1.25, linetype = 'dashed') +
      scale_colour_manual(values = rep('grey', max(plot_data$s))) +
      theme_bw() +
      theme(legend.position = 'none') +
      labs(x = input_name, y = output_name)
      
    if (length(inds) > 1){
      plot <- plot + facet_wrap(vars(.data$Run))
    }
    
    # If provided with true output to overlay
    if (!(is.null(Truth))){
      truth_data <- data.frame(Input = input_values, 
                               Output = c(Truth), 
                               s = 1,
                               Run = rep(inds, each = ell))
      plot <- plot +
        geom_line(data = truth_data, col = 'black', size = 1.25)
    }

  }
  
  else {
    ell <- dim(Samples$mean)[1]
    runs <- dim(Samples$mean)[2]
    
    # If only have a single input
    if (is.null(runs)){
      ell <- length(Samples$mean)
      inds <- 1
      runs <- 1
    }
    
    # If don't provide which inputs to plot, and there's more than 1 option, set as all
    if (is.null(inds)){
      inds <- 1:runs
    }
    
    # If have multiple inputs and a set of indices was provided
    if (runs > 1 & !(is.null(inds))){
      Samples$mean <- Samples$mean[,inds]
      Samples$lower <- Samples$lower[,inds]
      Samples$upper <- Samples$upper[,inds]
    }
    
    if (is.null(input_values)){
      input_values <- 1:ell
    }
    
    plot_data <- data.frame(Input = input_values, 
                            Output = c(Samples$mean),
                            Lower = c(Samples$lower),
                            Upper = c(Samples$upper),
                            Run = rep(inds, each = ell))
    
    plot <- ggplot(plot_data, aes(.data$Input, .data$Output)) + 
      geom_line(col = 'red', size = 1.25) +
      geom_line(aes(.data$Input, .data$Lower), col = 'red', size = 1.25, linetype = 'dashed') +
      geom_line(aes(.data$Input, .data$Upper), col = 'red', size = 1.25, linetype = 'dashed') +
      theme_bw() +
      theme(legend.position = 'none') +
      labs(x = input_name, y = output_name)
    
    if (length(inds) > 1){
      plot <- plot + facet_wrap(vars(.data$Run))
    }
    
    # If provided with true output to overlay
    if (!(is.null(Truth))){
      truth_data <- data.frame(Input = input_values, 
                               Output = c(Truth), 
                               Run = rep(inds, each = ell))
      plot <- plot +
        geom_line(data = truth_data, col = 'black', size = 1.25)
    }
  }
  return(plot)
}



#' Projects and reconstructs runs
#'
#' @param DataBasis 
#' @param q 
#' @param inds 
#' @param AddMean 
#' @param residual 
#' @param input_values 
#' @param input_name 
#' @param output_name 
#' @param ... 
#' 
#' @return
#' 
#' @export
PlotRecon <- function(DataBasis, q = 1, inds = 1:16, AddMean = TRUE, residual = FALSE, input_values = NULL, input_name = NULL, output_name = NULL, ...){
  
  ell <- nrow(DataBasis$tBasis)
  fields <- DataBasis$CentredField[,inds]
  basis <- DataBasis$tBasis[,1:q]
  k <- length(inds)
  recons <- lapply(1:k, function(i) ReconField(fields[,i], basis))
  recons <- matrix(unlist(recons), ell) 
  
  if (is.null(input_values)){
    input_values <- 1:ell
  }
  
  if (AddMean){
    mu <- DataBasis$EnsembleMean
  }
  
  else {
    mu <- 0*DataBasis$EnsembleMean
  }
  
  if (!(residual)){
    plot_data <- data.frame(Input = input_values, 
                            Output = c(recons) + mu,
                            Truth = c(fields) + mu,
                            Run = rep(inds, each = ell))
    
    plot <- ggplot(plot_data, aes(.data$Input, .data$Output)) + 
      geom_line(col = 'red', size = 1) +
      geom_line(aes(.data$Input, .data$Truth), col = 'black', size = 1) +
      facet_wrap(vars(.data$Run)) +
      theme_bw() +
      theme(legend.position = 'none') +
      labs(x = input_name, y = output_name)
  }
  
  if (residual){
    plot_data <- data.frame(Input = input_values, 
                            Residual = c(fields - recons),
                            Run = rep(inds, each = ell))
    
    plot <- ggplot(plot_data, aes(.data$Input, .data$Residual)) + 
      geom_line(col = 'red', size = 1) +
      facet_wrap(vars(.data$Run)) +
      theme_bw() +
      theme(legend.position = 'none') +
      labs(x = input_name, y = output_name)
  }

  return(plot)
}


#' Plots a 1D representation of the 1st q basis vectors
#'
#' @param DataBasis 
#' @param q 
#' @param input_values 
#' @param input_name 
#' 
#' @return 
#' 
#' @export
Plot1DBasis <- function(DataBasis, q = 9, input_values = NULL, input_name = NULL){
  ell <- nrow(DataBasis$tBasis)
  if (is.null(input_values)){
    input_values <- 1:ell
  }
  if (is.null(input_name)){
    input_name <- 'Input'
  }
  
  plot_basis <- data.frame(Values = input_values,
                           Vec = rep(1:q, each = ell),
                           Weight = c(DataBasis$tBasis[,1:q]))

  plot <- ggplot(plot_basis, aes(x = .data$Values, y = .data$Weight)) +
    geom_line() +
    facet_wrap(vars(.data$Vec)) +
    labs(x = input_name) +
    theme_bw() +
    theme(legend.position = 'none')
  
  return(plot)
}


#' Plots ensemble of 1D profile/TS
#'
#' @param DataBasis 
#' @param AddMean 
#' @param inds 
#' @param input_values 
#' @param input_name 
#' 
#' @return
#' 
#' @export
Plot1DData <- function(DataBasis, AddMean = TRUE, inds = NULL, input_values = NULL, input_name = NULL){
  ell <- nrow(DataBasis$tBasis)
  n <- ncol(DataBasis$CentredField)
  
  if (is.null(input_values)){
    input_values <- 1:ell
  }
  if (is.null(input_name)){
    input_name <- 'Input'
  }
  
  if (is.null(inds)){
    inds <- 1:n
  }
  
  if (AddMean){
    mu <- DataBasis$EnsembleMean
  }
  
  else {
    mu <- rep(0, ell)
  }
  
  plot_data <- data.frame(Values = input_values,
                          Output = c(DataBasis$CentredField[,inds] + mu),
                          Run = rep(inds, each = ell))
  
  plot <- ggplot(plot_data, aes(x = .data$Values, y = .data$Output, col = as.factor(.data$Run))) +
    geom_line() +
    labs(x = input_name) +
    theme_bw() +
    theme(legend.position = 'none')
  
  return(plot)
}



#' Plots proportion of variance cumulatively, or individually, explained by each basis vector
#'
#' @param DataBasis 
#' @param type 
#' @param ... 
#' 
#' @return 
#' 
#' @export
PlotExplained <- function(DataBasis, type = 'cumulative', ...){
  n <- ncol(DataBasis$tBasis)
  vars <- lapply(1:n, function(k) VarExplained(DataBasis$tBasis[,1:k], DataBasis$CentredField, ...))
  if (type == 'cumulative'){
    plot <- ggplot(data.frame(q = 1:n, Proportion = unlist(vars)), aes(x = .data$q, y = .data$Proportion)) +
      geom_line() +
      ylim(0,1) +
      labs(x = 'Vector') +
      theme_bw()
  }
  if (type == 'individual'){
    plot <- ggplot(data.frame(q = 1:n, Proportion = diff(c(0,unlist(vars)))), aes(x = .data$q, y = .data$Proportion)) +
      geom_point() +
      ylim(0,1) +
      labs(x = 'Vector') +
      theme_bw()
  }
  return(plot)
}


#' Plots reconstruction error for each basis truncation
#'
#' @param DataBasis 
#' @param obs 
#' @param qmax 
#' @param ... 
#' 
#' @return 
#' 
#' @export
PlotReconError <- function(DataBasis, obs, qmax = NULL, ...){
  n <- ncol(DataBasis$tBasis)
  if (is.null(qmax)){
    qmax <- n
  }
  RW <- errors(DataBasis$tBasis[,1:qmax], obs, ...)
  plot <- ggplot(data.frame(q = 1:qmax, y = RW), aes(x = .data$q, y = .data$y)) +
    geom_line(col = 'red') +
    geom_point(col = 'red', size = 0.75) +
    ylim(0,RW[1]) +
    labs(x = 'Vector', y = 'Reconstruction Error') +
    theme_bw() +
    geom_hline(yintercept = qchisq(0.995, length(obs))/length(obs), linetype = 'dashed')
  return(plot)
}


#' Plots residual for each basis truncation
#'
#' @param DataBasis 
#' @param obs 
#' @param q 
#' @param input_values 
#' @param input_name 
#' @param ... 
#'
#' @return
#' @export
PlotResid <- function(DataBasis, obs, q, input_values = NULL, input_name = NULL, ...){
  ell <- nrow(DataBasis$tBasis)
  if (is.null(input_values)){
    input_values <- 1:ell
  }
  if (is.null(input_name)){
    input_name <- 'Input'
  }
  
  plot_resids <- NULL
  for (i in 1:q){
    tmp <- ReconField(obs, DataBasis$tBasis[,1:i], ...)
    plot_resids <- rbind(plot_resids, data.frame(Input = input_values,
                                                 Error = obs - tmp,
                                                 q = i))
  }
  
  plot <- ggplot(plot_resids, aes(x = .data$Input, y = .data$Error)) +
    geom_line() +
    facet_wrap(vars(.data$q)) +
    labs(x = input_name) +
    theme_bw() +
    theme(legend.position = 'none')
  return(plot)
}






#' Plot pairs of coefficients and/or inputs, potentially coloured by a 3rd variable
#' 
#' Takes either output of `Project`, or something like tData (containing inputs and coefficients)
#'
#' @param coeffs 
#' @param x 
#' @param y 
#' @param col 
#' @param obs 
#' 
#' @return 
#' 
#' @export
PlotPair <- function(coeffs, x = 'C1', y = 'C2', col = NULL, obs = NULL){
  coeffs <- as.data.frame(coeffs)
  plot_data <- data.frame(x = coeffs[,x],
                          y = coeffs[,y])
  if (!(is.null(col))){
    plot_data$col <- coeffs[,col]
  }
  
  if (!(is.null(obs))){
    obs <- as.data.frame(obs)
    plot_data_obs <- data.frame(x = obs[,x],
                                y = obs[,y])
  }
  
  plot <- ggplot(plot_data, aes(.data$x, .data$y, col = .data$col)) +
    geom_point() +
    viridis::scale_colour_viridis() +
    labs(x = x, y = y, col = col) +
    theme_bw()
  
  if (!(is.null(obs))){
    plot <- plot + geom_point(data = plot_data_obs, col = 'red', size = 4, shape = 4)
  }

  return(plot)
}


#' Visualise active variables across set of emulators
#'
#' @param Ems 
#' @param InputNames 
#' 
#' @return 
#' 
#' @export
PlotActive <- function(Ems, InputNames){
  q <- length(Ems) # number of emulators
  plot_data <- NULL
  for (i in 1:q){
    plot_data <- rbind(plot_data, data.frame(Parameter = InputNames,
                                             Emulator = i,
                                             Active = InputNames %in% Ems[[i]]$active))
  }
  
  plot <- ggplot(plot_data, aes(x = as.factor(.data$Emulator), y = .data$Parameter, fill = as.factor(.data$Active))) +
    geom_tile() +
    labs(x = 'Emulator') +
    theme(legend.position = 'none') +
    scale_fill_manual(values = c('red', 'darkgreen'))
  
  return(plot)
}


#' Leave-one-out plotting in ggplot
#' 
#' @param emulator 
#' 
#' @return 
#' 
#' @export
LeaveOneOut <- function(emulator){
  em <- emulator$em
  loo_preds <- RobustGaSP::leave_one_out_rgasp(em)
  loo_preds$lower95 <- loo_preds$mean - 1.96*loo_preds$sd
  loo_preds$upper95 <- loo_preds$mean + 1.96*loo_preds$sd
  response <- emulator$train_data[,dim(emulator$train_data)[2]]
  loo_preds$truth <- response
  upp <- max(c(loo_preds$upper95, response))
  low <- min(c(loo_preds$lower95, response))
  
  loo_preds$In95 <- loo_preds$truth >= loo_preds$lower95 & loo_preds$truth <= loo_preds$upper95
  perc_outside <- round(sum(loo_preds$In95 == FALSE) / length(loo_preds$In95) * 100, 1)

  cols <- my.cols
  # Ensuring good points still coloured green if no points outside
  if (perc_outside == 0){
    cols[2:3] <- cols[3]
  }
  
  plot <- ggplot(as.data.frame(loo_preds), aes(x = .data$truth, y = .data$mean, col = .data$In95)) +
    geom_errorbar(aes(ymin = .data$lower95, ymax = .data$upper95), col = cols[1]) +
    geom_point() +
    scale_colour_manual(values = c(cols[2:3])) +
    geom_abline(slope = 1, alpha = 0.6) +
    labs(y = 'Prediction', x = 'Truth', title = paste0('Outside 95% = ', perc_outside, '%')) +
    theme_bw() +
    theme(legend.position = 'none')

  return(plot)
}


# Old version in base R plot
# LeaveOneOut <- function(emulator){
#   require(sfsmisc)
#   em <- emulator$em
#   loo_preds <- leave_one_out_rgasp(em)
#   loo_preds$lower95 <- loo_preds$mean - 1.96*loo_preds$sd
#   loo_preds$upper95 <- loo_preds$mean + 1.96*loo_preds$sd
#   response <- emulator$train_data[,dim(emulator$train_data)[2]]
#   upp <- max(c(loo_preds$upper95, response))
#   low <- min(c(loo_preds$lower95, response))
#   errbar(loo_preds$mean, loo_preds$mean, loo_preds$upper95, loo_preds$lower95, cap = 0.015, pch=20, 
#          ylim=c(low,upp),xlab = "Prediction",ylab="Data", main = 'Leave-one-out')
#   points(loo_preds$mean, response, pch=19,
#          col = ifelse(response > loo_preds$upper95 | response < loo_preds$lower95, "red", "green"))
# }




#' Validating a GaSP emulator
#' 
#' Given a validation dataset, predicts and plots the mean and 95% uncertainty interval against the true output
#' 
#' @param emulator either a single output from BuildGasp, or a list of emulators
#' @param ValidationData a validation data frame containing inputs and true output. If NULL, validation is performed using the validation_data output of the emulator
#' @param IndivPars 
#' 
#' @return 
#'  
#' @export
Validate <- function(emulator, ValidationData = NULL, IndivPars = FALSE){
  q <- length(emulator$em) # 1 if a single output from BuildGasp, 0 if a list of emulators
  if (q == 1){
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
    preds <- PredictGasp(design, emulator)
    upp <- max(c(preds$upper95, response))
    low <- min(c(preds$lower95, response))
    preds$truth <- response
    
    preds$In95 <- preds$truth >= preds$lower95 & preds$truth <= preds$upper95
    perc_outside <- round(sum(preds$In95 == FALSE) / length(preds$In95) * 100, 1)
    
    cols <- my.cols
    
    # Ensuring good points still coloured green if no points outside
    if (perc_outside == 0){
      cols[2:3] <- cols[3]
    }
    
    plots <- ggplot(as.data.frame(preds), aes(x = .data$truth, y = .data$mean, col = .data$In95)) +
      geom_errorbar(aes(ymin = .data$lower95, ymax = .data$upper95), col = cols[1]) +
      geom_point() +
      scale_colour_manual(values = c(cols[2:3])) +
      geom_abline(slope = 1, alpha = 0.6) +
      labs(y = 'Prediction', x = 'Truth', title = paste0('Outside 95% = ', perc_outside, '%')) +
      theme_bw() +
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
          theme_bw() +
          theme(legend.position = 'none')
      }
    }
  }
  else {
    
    plots <- NULL
    plots_inputs <- NULL
    
    for (i in 1:length(emulator)){
      if (!is.null(ValidationData)){
        emulator[[i]]$validation_data <- ValidationData[,colnames(emulator[[i]]$train_data)]
      }
      resp_ind <- dim(emulator[[i]]$validation_data)[2]
      design <- emulator[[i]]$validation_data[,-resp_ind]
      if (length(emulator[[i]]$active) == 1){
        design <- as.matrix(design, ncol = 1)
        colnames(design) <- emulator[[i]]$active
      }
      response <- emulator[[i]]$validation_data[,resp_ind]
      preds <- PredictGasp(design, emulator[[i]])
      upp <- max(c(preds$upper95, response))
      low <- min(c(preds$lower95, response))
      preds$truth <- response
      
      preds$In95 <- preds$truth >= preds$lower95 & preds$truth <= preds$upper95
      perc_outside <- round(sum(preds$In95 == FALSE) / length(preds$In95) * 100, 1)
      cols <- my.cols
      if (perc_outside == 0){
        cols[2:3] <- cols[3]
      }
      
      plots[[i]] <- ggplot(as.data.frame(preds), aes(x = .data$truth, y = .data$mean, col = .data$In95)) +
        geom_errorbar(aes(ymin = .data$lower95, ymax = .data$upper95), col = cols[1]) +
        geom_point() +
        scale_colour_manual(values = c(cols[2:3])) +
        geom_abline(slope = 1, alpha = 0.6) +
        labs(y = 'Prediction', x = 'Truth', title = paste0('Outside 95% = ', perc_outside, '%')) +
        theme_bw() +
        theme(legend.position = 'none')

      if (IndivPars == TRUE){
        plots_inputs[[i]] <- list()
        for (j in 1:dim(design)[2]){
          plot_data <- data.frame(x = design[,j], as.data.frame(preds))
          plots_inputs[[i]][[j]] <- ggplot(plot_data, aes(x = .data$x, y = .data$mean, col = .data$In95)) +
            geom_errorbar(aes(ymin = .data$lower95, ymax = .data$upper95), col = cols[1]) +
            geom_point() +
            scale_colour_manual(values = c(cols[2:3])) +
            labs(y = 'Prediction', x = paste0(colnames(design)[j])) +
            theme_bw() +
            theme(legend.position = 'none')
        }
      }
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

ValidateGasp <- Validate

#' Old version in base R
# ValidateGasp <- function(emulator, ValidationData = NULL, IndivPars = FALSE){
#   require(sfsmisc)
#   q <- length(emulator$em) # 1 if a single output from BuildGasp, 0 if a list of emulators
#   if (q == 1){
#     if (!is.null(ValidationData)){
#       emulator$validation_data <- ValidationData[,colnames(emulator$train_data)]
#     }
#     resp_ind <- dim(emulator$validation_data)[2]
#     design <- emulator$validation_data[,-resp_ind]
#     if (length(emulator$active) == 1){
#       design <- as.matrix(design, ncol = 1)
#       colnames(design) <- emulator$active
#     }
#     response <- emulator$validation_data[,resp_ind]
#     preds <- PredictGasp(design, emulator)
#     upp <- max(c(preds$upper95, response))
#     low <- min(c(preds$lower95, response))
#     errbar(preds$mean, preds$mean, preds$upper95, preds$lower95, cap = 0.015, pch=20, 
#            ylim=c(low,upp), main="",xlab = "Prediction",ylab="Data")
#     points(preds$mean, response, pch=19,
#            col = ifelse(response > preds$upper95 | response < preds$lower95, "red", "green"))
#     if (IndivPars == TRUE){
#       for (i in 1:dim(design)[2]){
#         errbar(design[,i], preds$mean, preds$upper95, preds$lower95, cap = 0.015, pch=20, 
#                ylim=c(low,upp), xlab = "Input",ylab="Prediction", main = paste(colnames(design)[i]))
#         points(design[,i], response, pch=19,
#                col = ifelse(response > preds$upper95 | response < preds$lower95, "red", "green"))
#       }
#     }
#   }
#   else {
#     for (i in 1:length(emulator)){
#       if (!is.null(ValidationData)){
#         emulator[[i]]$validation_data <- ValidationData[,colnames(emulator[[i]]$train_data)]
#       }
#       resp_ind <- dim(emulator[[i]]$validation_data)[2]
#       design <- emulator[[i]]$validation_data[,-resp_ind]
#       if (length(emulator[[i]]$active) == 1){
#         design <- as.matrix(design, ncol = 1)
#         colnames(design) <- emulator[[i]]$active
#       }
#       response <- emulator[[i]]$validation_data[,resp_ind]
#       preds <- PredictGasp(design, emulator[[i]])
#       upp <- max(c(preds$upper95, response))
#       low <- min(c(preds$lower95, response))
#       errbar(preds$mean, preds$mean, preds$upper95, preds$lower95, cap = 0.015, pch=20, 
#              ylim=c(low,upp),xlab = "Prediction",ylab="Data", main = colnames(emulator[[i]]$validation_data)[resp_ind])
#       points(preds$mean, response, pch=19,
#              col = ifelse(response > preds$upper95 | response < preds$lower95, "red", "green"))
#       if (IndivPars == TRUE){
#         for (j in 1:dim(design)[2]){
#           errbar(design[,j], preds$mean, preds$upper95, preds$lower95, cap = 0.015, pch=20, 
#                  ylim=c(low,upp),xlab = "Input",ylab="Prediction", main = paste(colnames(design)[j]))
#           points(design[,j], response, pch=19,
#                  col = ifelse(response > preds$upper95 | response < preds$lower95, "red", "green"))
#         }
#       }
#     }
#   }
# }





#' Pairs plot of NROY space
#'
#' @param Design Data frame of inputs, possibly with a TRUE/FALSE final column named `NROY`.
#' @param k Indices relating to which columns of `Design` to plot. If NULL, plots all.
#' @param NROY Vector of TRUE/FALSE labels, corresponding to classification of `Design`. If NULL, uses the final column of `Design`.
#' @param Truth NOT YET ADDED
#' @param size Controlling size of points on lower half of plot
#'
#' @return
#' @export
#'
#' @examples
PlotNROY <- function(Design, k = NULL, NROY = NULL, Truth = NULL, size = 0.5){
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
  
  plot <- GGally::ggpairs(Design, columns = k, aes(col = NROY), 
                          upper = list(continuous = GGally::wrap("density", alpha = 0.5), combo = "box_no_facet"),
                          lower = list(continuous = GGally::wrap("points", alpha = 0.3, size = size), combo = GGally::wrap("dot_no_facet", alpha = 0.4)),
                          diag = list(continuous = GGally::wrap("densityDiag", alpha = 0.3)), 
                          legend = 1) +
    theme_bw() +
    theme(legend.position = "bottom") + 
    scale_colour_manual(values = my.cols[-1]) +
    scale_fill_manual(values = my.cols[-1])
  
  return(plot)
}



