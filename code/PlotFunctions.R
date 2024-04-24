# Some additional basis functionality, that will be combined with other source code at some point
# Mostly plotting, but also sampling from emulator posteriors

library(ggplot2)
library(reshape2)
library(cowplot)

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


#' Takes samples (or a summary of these) from BasisEmSamples and plots, adds the truth if this is provided
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
    
    plot <- ggplot(plot_data, aes(Input, Output, col = as.factor(s))) + 
      geom_line() +
      geom_line(data = data.frame(aggregate(Output ~ Input + Run, plot_data, mean)), col = 'red', size = 1.25) +
      geom_line(data = data.frame(aggregate(Output ~ Input + Run, plot_data, quantile, probs = 0.025)), col = 'red', size = 1.25, linetype = 'dashed') +
      geom_line(data = data.frame(aggregate(Output ~ Input + Run, plot_data, quantile, probs = 0.975)), col = 'red', size = 1.25, linetype = 'dashed') +
      scale_colour_manual(values = rep('grey', max(plot_data$s))) +
      theme(legend.position = 'none') +
      labs(x = input_name, y = output_name)
      
    if (length(inds) > 1){
      plot <- plot + facet_wrap(vars(Run))
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
    
    plot <- ggplot(plot_data, aes(Input, Output)) + 
      geom_line(col = 'red', size = 1.25) +
      geom_line(aes(Input, Lower), col = 'red', size = 1.25, linetype = 'dashed') +
      geom_line(aes(Input, Upper), col = 'red', size = 1.25, linetype = 'dashed') +
      theme(legend.position = 'none') +
      labs(x = input_name, y = output_name)
    
    if (length(inds) > 1){
      plot <- plot + facet_wrap(vars(Run))
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
    
    plot <- ggplot(plot_data, aes(Input, Output)) + 
      geom_line(col = 'red', size = 1) +
      geom_line(aes(Input, Truth), col = 'black', size = 1) +
      facet_wrap(vars(Run)) +
      theme(legend.position = 'none') +
      labs(x = input_name, y = output_name)
  }
  
  if (residual){
    plot_data <- data.frame(Input = input_values, 
                            Residual = c(fields - recons),
                            Run = rep(inds, each = ell))
    
    plot <- ggplot(plot_data, aes(Input, Residual)) + 
      geom_line(col = 'red', size = 1) +
      facet_wrap(vars(Run)) +
      theme(legend.position = 'none') +
      labs(x = input_name, y = output_name)
  }

  return(plot)
}


#' Plots a 1D representation of the 1st q basis vectors
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

  plot <- ggplot(plot_basis, aes(x = Values, y = Weight)) +
    geom_line() +
    facet_wrap(vars(Vec)) +
    labs(x = input_name) +
    theme(legend.position = 'none')
  
  return(plot)
}

#' Plots proportion of variance cumulatively, or individually, explained by each basis vector
PlotExplained <- function(DataBasis, type = 'cumulative', ...){
  n <- ncol(DataBasis$tBasis)
  vars <- lapply(1:n, function(k) VarExplained(DataBasis$tBasis[,1:k], DataBasis$CentredField, ...))
  if (type == 'cumulative'){
    plot <- ggplot(data.frame(q = 1:n, Proportion = unlist(vars)), aes(x = q, y = Proportion)) +
      geom_line() +
      ylim(0,1) +
      labs(x = 'Vector')
  }
  if (type == 'individual'){
    plot <- ggplot(data.frame(q = 1:n, Proportion = diff(c(0,unlist(vars)))), aes(x = q, y = Proportion)) +
      geom_point() +
      ylim(0,1) +
      labs(x = 'Vector')
  }
  return(plot)
}


#' Plots reconstruction error for each basis truncation
PlotReconError <- function(DataBasis, obs, qmax = NULL, ...){
  n <- ncol(DataBasis$tBasis)
  if (is.null(qmax)){
    qmax <- n
  }
  RW <- errors(DataBasis$tBasis[,1:qmax], obs, ...)
  plot <- ggplot(data.frame(q = 1:qmax, y = RW), aes(x = q, y = y)) +
    geom_line(col = 'red') +
    geom_point(col = 'red', size = 0.75) +
    ylim(0,RW[1]) +
    labs(x = 'Vector', y = 'Reconstruction Error') +
    geom_hline(yintercept = qchisq(0.995, length(obs))/length(obs), linetype = 'dashed')
  return(plot)
}


#' Plot pairs of coefficients and/or inputs, potentially coloured by a 3rd variable
#' Takes either output of `Project`, or something like tData (containing inputs and coefficients)
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
  
  plot <- ggplot(plot_data, aes(x, y, col = col)) +
    geom_point() +
    scale_colour_viridis() +
    labs(x = x, y = y, col = col)
  
  if (!(is.null(obs))){
    plot <- plot + geom_point(data = plot_data_obs, col = 'red', size = 4, shape = 4)
  }

  return(plot)
}


#' Visualise active variables across set of emulators
PlotActive <- function(Ems, InputNames){
  q <- length(Ems) # number of emulators
  plot_data <- NULL
  for (i in 1:q){
    plot_data <- rbind(plot_data, data.frame(Parameter = InputNames,
                                             Emulator = i,
                                             Active = InputNames %in% Ems[[i]]$active))
  }
  
  plot <- ggplot(plot_data, aes(x = as.factor(Emulator), y = Parameter, fill = as.factor(Active))) +
    geom_tile() +
    labs(x = 'Emulator') +
    theme(legend.position = 'none') +
    scale_fill_manual(values = c('red', 'darkgreen'))
  
  return(plot)
}


#' Leave-one-out plotting in ggplot
#' 
LeaveOneOut <- function(emulator){
  em <- emulator$em
  loo_preds <- leave_one_out_rgasp(em)
  loo_preds$lower95 <- loo_preds$mean - 1.96*loo_preds$sd
  loo_preds$upper95 <- loo_preds$mean + 1.96*loo_preds$sd
  response <- emulator$train_data[,dim(emulator$train_data)[2]]
  loo_preds$truth <- response
  upp <- max(c(loo_preds$upper95, response))
  low <- min(c(loo_preds$lower95, response))
  
  loo_preds$In95 <- loo_preds$truth >= loo_preds$lower95 & loo_preds$truth <= loo_preds$upper95
  perc_outside <- round(sum(loo_preds$In95 == FALSE) / length(loo_preds$In95) * 100, 1)
  cols <- c('darkgrey', viridis(100)[31], viridis(100)[81])
  
  plot <- ggplot(as.data.frame(loo_preds), aes(x = truth, y = mean, col = In95)) +
    geom_errorbar(aes(ymin = lower95, ymax = upper95), col = cols[1]) +
    geom_point() +
    scale_colour_manual(values = c(cols[2:3])) +
    geom_abline(slope = 1, alpha = 0.6) +
    labs(y = 'Prediction', x = 'Truth', title = paste0('Outside 95% = ', perc_outside, '%')) +
    theme(legend.position = 'none')
  
  return(plot)
}


#' Prediction over validation set in ggplot
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
    cols <- c('darkgrey', viridis(100)[31], viridis(100)[81])
    
    plots <- ggplot(as.data.frame(preds), aes(x = truth, y = mean, col = In95)) +
      geom_errorbar(aes(ymin = lower95, ymax = upper95), col = cols[1]) +
      geom_point() +
      scale_colour_manual(values = c(cols[2:3])) +
      geom_abline(slope = 1, alpha = 0.6) +
      labs(y = 'Prediction', x = 'Truth', title = paste0('Outside 95% = ', perc_outside, '%')) +
      theme(legend.position = 'none')

    if (IndivPars == TRUE){
      plots_inputs <- NULL
      for (j in 1:dim(design)[2]){
        plot_data <- data.frame(x = design[,j], as.data.frame(preds))
        plots_inputs[[j]] <- ggplot(plot_data, aes(x = x, y = mean, col = In95)) +
          geom_errorbar(aes(ymin = lower95, ymax = upper95), col = cols[1]) +
          geom_point() +
          scale_colour_manual(values = c(cols[2:3])) +
          labs(y = 'Prediction', x = paste0(colnames(design)[j])) +
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
      cols <- c('darkgrey', viridis(100)[31], viridis(100)[81])
      
      plots[[i]] <- ggplot(as.data.frame(preds), aes(x = truth, y = mean, col = In95)) +
        geom_errorbar(aes(ymin = lower95, ymax = upper95), col = cols[1]) +
        geom_point() +
        scale_colour_manual(values = c(cols[2:3])) +
        geom_abline(slope = 1, alpha = 0.6) +
        labs(y = 'Prediction', x = 'Truth', title = paste0('Outside 95% = ', perc_outside, '%')) +
        theme(legend.position = 'none')

      if (IndivPars == TRUE){
        plots_inputs[[i]] <- list()
        for (j in 1:dim(design)[2]){
          plot_data <- data.frame(x = design[,j], as.data.frame(preds))
          plots_inputs[[i]][[j]] <- ggplot(plot_data, aes(x = x, y = mean, col = In95)) +
            geom_errorbar(aes(ymin = lower95, ymax = upper95), col = cols[1]) +
            geom_point() +
            scale_colour_manual(values = c(cols[2:3])) +
            labs(y = 'Prediction', x = paste0(colnames(design)[j])) +
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





