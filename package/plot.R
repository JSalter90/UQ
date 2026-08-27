#' Functions for plotting

#' Plots the first q basis vectors
#' 
#' 
#'
#' @param DataBasis 
#' @param q 
#' @param output_inds
#' 
#' @returns description
#' 
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_errorbar scale_colour_manual labs vars geom_abline theme theme_bw facet_wrap
#' @importFrom rlang .data
NULL
#' 
#' @export
PlotBasis <- function(DataBasis, q = 9, output_inds = NULL){
  
  # Dimension of basis
  ell <- nrow(DataBasis$tBasis)
  
  # If a specific subset hasn't been selected, use all ell dimensions
  if (is.null(output_inds)){
    output_inds <- 1:ell
  }
  
  if (any(output_inds > ell)){
    stop('Selected output_inds contain indices larger than the size of the data')
  }
  
  # Extract dimension names and locations from DataBasis
  input_names <- colnames(DataBasis$grid)
  input_values <- DataBasis$grid[output_inds,]
  
  # If not provided, set
  if (is.null(input_names)){
    input_names <- 'Input'
  }
  
  if (is.null(input_values)){
    input_values <- output_inds
  }

  # Flag for whether the 2nd dimension relates to grid values or multiple variables
  if (length(input_names) > 1){
    multi_var <- ifelse(length(unique(DataBasis$grid[,2])) < 10, TRUE, FALSE)
  }

  if (length(input_names) == 1){
    plot_data <- data.frame(Input = input_values,
                            Vector = rep(1:q, each = length(output_inds)),
                            Weight = c(DataBasis$tBasis[output_inds,1:q]))
    
    plot <- ggplot(plot_data, aes(x = .data$Input, y = .data$Weight)) +
      geom_line() +
      facet_wrap(vars(.data$Vector)) +
      labs(x = input_names) +
      theme_bw() +
      theme(legend.position = 'none')
  }
  
  else if (length(input_names) == 2 & !multi_var){
    plot_data <- data.frame(Input1 = input_values[,1],
                            Input2 = input_values[,2],
                            Vector = rep(1:q, each = length(output_inds)),
                            Weight = c(DataBasis$tBasis[output_inds,1:q]))
    
    point_size <- ifelse(length(output_inds) >= 1000 & q > 3, 0.25, 0.75)
    
    plot <- ggplot(plot_data, aes(.data$Input1, .data$Input2, col = .data$Weight)) +
      geom_point(size = point_size) +
      facet_wrap(vars(.data$Vector)) +
      labs(x = input_names[1], y = input_names[2]) +
      viridis::scale_colour_viridis() +
      theme_bw()
  }
  
  else if (length(input_names) == 2 & multi_var){
    plot_data <- data.frame(Input = input_values[,1],
                            Variable = input_values[,2],
                            Vector = rep(1:q, each = length(output_inds)),
                            Weight = c(DataBasis$tBasis[output_inds,1:q]))
    
    plot <- ggplot(plot_data, aes(.data$Input, .data$Weight)) +
      geom_line() +
      facet_grid(.data$Variable ~ .data$Vector) +
      labs(x = input_names[1]) +
      theme_bw()
  }
  
  else if (length(input_names) == 3){
    
    
  }

  return(plot)
}




# Extend to 2D, in which case need to select inds

# Adding obs only works where 1D profiles, not 2D output. But could add as a panel?

#' Plots ensemble of 1D profile/TS
#'
#' @param DataBasis 
#' @param AddMean 
#' @param data_inds Which (simulation) data to plot. If NULL, plots all contained in `DataBasis`.
#' @param output_inds Which output locations to plot. If NULL, plots all.
#' 
#' @returns description
#' 
#' @export
PlotData <- function(DataBasis, AddMean = TRUE, data_inds = NULL, output_inds = NULL, Obs = NULL, MeanOnly = FALSE){
  
  ell <- nrow(DataBasis$tBasis)
  n <- ncol(DataBasis$CentredField)
  
  # If a specific subset hasn't been selected, use all ell dimensions
  if (is.null(output_inds)){
    output_inds <- 1:ell
  }
  
  if (any(output_inds > ell)){
    stop('Selected output_inds contain indices larger than the size of the data')
  }
  
  # If a specific subset hasn't been selected, plot all n outputs in DataBasis$CentredField
  if (is.null(data_inds)){
    data_inds <- 1:n
  }
  
  if (any(data_inds > n)){
    stop('Selected data_inds contain indices larger than the number of available simulations')
  }
  
  # Extract dimension names and locations from DataBasis
  input_values <- DataBasis$grid[output_inds,]
  input_names <- colnames(DataBasis$grid)
  
  # If not provided, set
  if (is.null(input_names)){
    input_names <- 'Input'
  }
  
  if (is.null(input_values)){
    input_values <- output_inds
  }

  if (MeanOnly){
    DataBasis$CentredField <- 0*DataBasis$CentredField
    mu <- DataBasis$EnsembleMean
    data_inds <- 1
  }
  
  else if (AddMean){
    mu <- DataBasis$EnsembleMean[output_inds]
  }
  
  else {
    mu <- rep(0, length(output_inds))
  }
  
  # Process obs if provided
  if (!(is.null(Obs))){
    Obs <- Obs[output_inds]
    
    if (!(AddMean)){
      Obs <- Obs - DataBasis$EnsembleMean[output_inds]
    }
  }

  # Flag for whether the 2nd dimension relates to grid values or multiple variables
  if (length(input_names) > 1){
    multi_var <- ifelse(length(unique(DataBasis$grid[,2])) < 10, TRUE, FALSE)
  }

  if (length(input_names) == 1){
    plot_data <- data.frame(Input = input_values,
                            Output = c(DataBasis$CentredField[output_inds,data_inds] + mu),
                            Run = rep(data_inds, each = length(output_inds)))
    
    if (MeanOnly){
      plot_data$Run <- 'Mean'
    }
    
    plot <- ggplot(plot_data, aes(x = .data$Input, y = .data$Output, col = as.factor(.data$Run))) +
      geom_line() +
      labs(x = input_names) +
      theme_bw() +
      theme(legend.position = 'none')
    
    if (!(is.null(Obs))){
      Obs <- data.frame(Input = input_values,
                        Output = Obs)
      plot <- plot + geom_line(data = Obs, aes(x = .data$Input, y = .data$Output), col = cols_validate[2], size = 1.5)
    }
    
  }
  
  else if (length(input_names) == 2 & !multi_var){
    plot_data <- data.frame(Input1 = input_values[,1],
                            Input2 = input_values[,2],
                            Output = c(DataBasis$CentredField[output_inds,data_inds] + mu[output_inds]),
                            Run = rep(data_inds, each = length(output_inds)))
    
    if (MeanOnly){
      plot_data$Run <- 'Mean'
    }
    
    point_size <- ifelse(length(output_inds) >= 1000 & length(data_inds) > 3, 0.25, 0.75)
    
    plot <- ggplot(plot_data, aes(x = .data$Input1, y = .data$Input2, col = .data$Output)) +
      geom_point(size = point_size) +
      labs(x = input_names) +
      facet_wrap(vars(.data$Run), scales = 'fixed') +
      viridis::scale_colour_viridis() +
      theme_bw()
  }
  
  else if (length(input_names) == 2 & multi_var){
    plot_data <- data.frame(Input = input_values[,1],
                            Variable = input_values[,2],
                            Output = c(DataBasis$CentredField[output_inds,data_inds] + mu[output_inds]),
                            Run = rep(data_inds, each = length(output_inds)))
    
    if (MeanOnly){
      plot_data$Run <- 'Mean'
    }
    
    plot <- ggplot(plot_data, aes(x = .data$Input, y = .data$Output, col = as.factor(.data$Run))) +
      geom_line() +
      labs(x = input_names) +
      facet_wrap(vars(.data$Variable)) +
      theme_bw() +
      theme(legend.position = 'none')
    
    if (!(is.null(Obs))){
      Obs <- data.frame(Input = input_values[,1],
                        Variable = input_values[,2],
                        Output = Obs)
      plot <- plot + geom_line(data = Obs, aes(x = .data$Input, y = .data$Output), col = cols_validate[2], size = 1.5)
    }
  }
  
  else if (length(input_names) == 3){
    
  }
  
  
  
  return(plot)
}




# Currently just does truncations for obs, but surely could do any ensemble member
# Similar idea to PlotData
# Could tag as a 'variable' the number used to recon
# Could plot full ensemble of residuals given q

#' Plots residual for each basis truncation
#'
#' @param DataBasis 
#' @param obs 
#' @param q 
#' @param input_values 
#' @param input_name 
#' @param ... 
#'
#' @returns 
#' 
#' @export
PlotResid <- function(DataBasis, obs, q, ...){
  
  ell <- nrow(DataBasis$tBasis)
  #n <- ncol(DataBasis$CentredField)
  
  # Extract dimension names and locations from DataBasis
  input_values <- DataBasis$grid[1:ell,]
  input_names <- colnames(DataBasis$grid)
  
  # If not provided, set
  if (is.null(input_names)){
    input_names <- 'Input'
  }
  
  if (is.null(input_values)){
    input_values <- 1:ell
  }
  
  # Option 1 - plot all residuals for given q
  
  
  
  
  # Option 2 - plot all residuals for multiple q
  
  # Option 3 - plot residual for single point
  
  
  
  
  
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
    labs(x = input_names[1]) +
    theme_bw() +
    theme(legend.position = 'none')
  return(plot)
}



# similar plotting to PlotData?
# need to allow alternative fields?

#' Projects and reconstructs runs
#'
#' @param DataBasis 
#' @param q 
#' @param data_inds 
#' @param AddMean 
#' @param residual 
#' @param input_values 
#' @param input_name 
#' @param output_name 
#' @param ... 
#' 
#' @returns description
#' 
#' @export
PlotRecon <- function(DataBasis, q = 1, data_inds = 1:16, AddMean = TRUE, residual = FALSE, input_values = NULL, input_name = NULL, output_name = NULL, ...){
  
  ell <- nrow(DataBasis$tBasis)
  n <- ncol(DataBasis$CentredField)
  
  if (any(data_inds > n)){
    stop('Selected data_inds contain indices larger than the number of available simulations')
  }
  
  fields <- DataBasis$CentredField[,data_inds]
  basis <- DataBasis$tBasis[,1:q]
  recons <- ReconField(fields, basis)
  
  # Extract dimension names and locations from DataBasis
  input_values <- DataBasis$grid[1:ell,]
  input_names <- colnames(DataBasis$grid)
  
  # If not provided, set
  if (is.null(input_names)){
    input_names <- 'Input'
  }
  
  if (is.null(input_values)){
    input_values <- 1:ell
  }
  
  if (is.null(output_name)){
    output_name <- 'Output'
  }
  
  if (AddMean){
    mu <- DataBasis$EnsembleMean
  }
  
  else {
    mu <- 0*DataBasis$EnsembleMean
  }
  
  # Flag for whether the 2nd dimension relates to grid values or multiple variables
  if (length(input_names) > 1){
    multi_var <- ifelse(length(unique(DataBasis$grid[,2])) < 10, TRUE, FALSE)
  }
  
  if (length(input_names) == 1){
    if (!(residual)){
      plot_data <- data.frame(Input = input_values, 
                              Output = c(c(recons) + mu, c(fields) + mu),
                              Type = rep(c('Recon', 'Truth'), each = ell*length(data_inds)),
                              Run = rep(data_inds, each = ell))
      
      plot <- ggplot(plot_data, aes(.data$Input, .data$Output, col = .data$Type)) + 
        geom_line(linewidth = 1) +
        scale_colour_manual(values = c('red', 'black')) +
        facet_wrap(vars(.data$Run)) +
        theme_bw() +
        labs(x = input_names[1], y = output_name)
    }
    
    if (residual){
      plot_data <- data.frame(Input = input_values, 
                              Residual = c(fields - recons),
                              Run = rep(data_inds, each = ell))
      
      plot <- ggplot(plot_data, aes(.data$Input, .data$Residual)) + 
        geom_line(col = 'red', linewidth = 1) +
        facet_wrap(vars(.data$Run)) +
        theme_bw() +
        theme(legend.position = 'none') +
        labs(x = input_names[1], y = output_name)
    }
  }
  
  else if (length(input_names) == 2 & !multi_var){
    if (!(residual)){
      plot_data <- data.frame(Input1 = input_values[,1],
                              Input2 = input_values[,2],
                              Output = c(recons) + mu,
                              Run = rep(data_inds, each = ell))
      
      point_size <- ifelse(ell >= 1000 & length(data_inds) > 3, 0.25, 0.75)
      
      plot <- ggplot(plot_data, aes(x = .data$Input1, y = .data$Input2, col = .data$Output)) + 
        geom_point(size = point_size) +
        viridis::scale_colour_viridis() +
        facet_wrap(vars(.data$Run), scales = 'fixed') +
        theme_bw() +
        labs(x = input_names[1], y = input_names[2])
    }
    
    if (residual){
      plot_data <- data.frame(Input1 = input_values[,1],
                              Input2 = input_values[,2],
                              Residual = c(fields - recons),
                              Run = rep(data_inds, each = ell))
      
      point_size <- ifelse(ell >= 1000 & length(data_inds) > 3, 0.25, 0.75)
      
      plot <- ggplot(plot_data, aes(x = .data$Input1, y = .data$Input2, col = .data$Residual)) + 
        geom_point(size = point_size) +
        ggplot2::scale_colour_gradient2(low = 'blue', mid = 'white', high = 'red') +
        facet_wrap(vars(.data$Run), scales = 'fixed') +
        theme_bw() +
        labs(x = input_names[1], y = input_names[2])
    }
    
    
    
    
  }
  
  else if (length(input_names) == 2 & multi_var){
    if (!(residual)){
      plot_data <- data.frame(Input = input_values[,1],
                              Variable = input_values[,2],
                              Output = c(c(recons) + mu, c(fields) + mu),
                              Type = rep(c('Recon', 'Truth'), each = ell*length(data_inds)),
                              Run = rep(data_inds, each = ell))
      
      plot <- ggplot(plot_data, aes(.data$Input, .data$Output, col = .data$Type)) + 
        geom_line(linewidth = 1) +
        scale_colour_manual(values = c('red', 'black')) +
        facet_grid(.data$Variable ~ .data$Run) +
        theme_bw() +
        labs(x = input_names[1], y = output_name)
    }
    
    if (residual){
      plot_data <- data.frame(Input = input_values[,1],
                              Variable = input_values[,2],
                              Residual = c(fields - recons),
                              Run = rep(data_inds, each = ell))
      
      plot <- ggplot(plot_data, aes(.data$Input, .data$Residual)) + 
        geom_line(col = 'red', linewidth = 1) +
        facet_grid(.data$Variable ~ .data$Run) +
        theme_bw() +
        theme(legend.position = 'none') +
        labs(x = input_names[1], y = output_name)
    }
  }
  
  return(plot)
}








#' Plotting emulator samples
#'
#' Takes samples (or a summary of these) from BasisEmSamples and plots, adds the truth if this is provided
#'
#' @param output_name 
#' @param sample_inds Only plot a subset of the provided samples. Still uses all samples to calculate mean/95%.
#' @param DataBasis 
#' @param emulator 
#' @param design 
#' @param samples 
#' @param data_inds 
#' @param output_inds 
#' @param interval 
#' @param AddMean 
#' @param Obs 
#' @param ... 
#' 
#' @returns
#' 
#' @export
PlotSamples <- function(DataBasis, emulator = NULL, design = NULL, samples = NULL, output_name = NULL, data_inds = NULL, sample_inds = NULL, output_inds = NULL, interval = 0.95, AddMean = TRUE, Obs = NULL,...){

  # Various combinations of design, samples, emulator can be provided
  if (!(is.null(design)) & is.null(samples) & !(is.null(data_inds))){  
    # Generate emulator predictions at required locations only
    preds <- Predict(emulator, design[data_inds,])
    # Draw samples
    Samples <- BasisEmSamples(preds, DataBasis, AddMean = AddMean, ...)
  }
  
  else if (is.null(design) & is.null(samples)){
    stop('Must provide either a design and emulators for prediction, or a set of samples')
  }
  
  # If data_inds not provided, use all of design
  else if (is.null(samples) & is.null(data_inds)){
    preds <- Predict(emulator, design)
    Samples <- BasisEmSamples(preds, DataBasis, AddMean = AddMean, ...)
  }
  
  else if (is.null(design) & !(is.null(samples))){
    Samples <- samples
  }
  
  if (interval <= 0 | interval >= 1){
    stop('interval must be between 0 and 1')
  }
  
  interval <- c((1-interval)/2, 1 - (1-interval)/2)

  # Construct truth
  if (!(is.null(data_inds))){
    if (AddMean){
      Truth <- DataBasis$CentredField[,data_inds] + DataBasis$EnsembleMean
    }
    else {
      Truth <- DataBasis$CentredField[,data_inds]
    }
    
    inds <- data_inds
  }
  
  else {
    tmp <- ifelse(is.na(dim(Samples)[3]), 1, dim(Samples)[3])
    inds <- 1:tmp
    Truth <- NULL
  }

  ell <- nrow(DataBasis$tBasis) # Dimension of field
  ns <- dim(Samples)[2] # number of samples

  # If don't provide a subset to plot, use all
  if (is.null(output_inds)){
    output_inds <- 1:ell
  }
  
  if (any(output_inds > ell)){
    stop('Selected output_inds contain indices larger than the size of the data')
  }
  
  # Extract dimension names and locations from DataBasis
  input_values <- DataBasis$grid[output_inds,]
  input_names <- colnames(DataBasis$grid)
  
  # If not provided, set
  if (is.null(input_names)){
    input_names <- 'Input'
  }
  
  if (is.null(input_values)){
    input_values <- output_inds
  }
  
  if (is.null(output_name)){
    output_name <- 'Output'
  }
  
  # Number of input locations
  n <- length(inds)

  # If have multiple inputs and a set of locations was provided
  if (length(inds) > 1 & !(is.null(output_inds))){
    Samples <- Samples[output_inds,,]
    if (!(is.null(Truth))){
      Truth <- Truth[output_inds,]
    }
  }
  
  # If have a single input
  if (length(inds) == 1 & !(is.null(output_inds))){
    Samples <- Samples[output_inds,]
    if (!(is.null(Truth))){
      Truth <- Truth[output_inds]
    }
  }

  # If didn't provide subset of samples to plot in background, plot all
  if (is.null(sample_inds)){
    sample_inds <- 1:ns
  }
  
  # If Obs provided, select plotting inds
  if (!(is.null(Obs))){
    Obs <- Obs[output_inds]
  }
  
  
  # Flag for whether the 2nd dimension relates to grid values or multiple variables
  if (length(input_names) > 1){
    multi_var <- ifelse(length(unique(DataBasis$grid[,2])) < 10, TRUE, FALSE)
  }
  
  if (length(input_names) == 1){
    plot_data <- data.frame(Input = input_values,
                            Output = c(Samples), 
                            s = rep(1:ns, each = length(output_inds)),
                            Run = rep(inds, each = length(output_inds)*ns),
                            Type = 'Samples')
    
    # Calculate summaries using all samples
    plot_data_summary <- rbind(data.frame(stats::aggregate(Output ~ Input + Run, plot_data, mean), Type = 'Mean'),
                               data.frame(stats::aggregate(Output ~ Input + Run, plot_data, stats::quantile, probs = interval[1]), Type = 'Lower'),
                               data.frame(stats::aggregate(Output ~ Input + Run, plot_data, stats::quantile, probs = interval[2]), Type = 'Upper'))
    
    plot_data_summary$s <- max(plot_data$s)+1
    plot_data_summary$s[which(plot_data_summary$Type == 'Lower')] <- max(plot_data$s)+2 # to allow to set different linetype
    plot_data_summary$s[which(plot_data_summary$Type == 'Upper')] <- max(plot_data$s)+3 # to allow to set different linetype
    plot_data_summary$Type[which(plot_data_summary$Type %in% c('Lower','Upper'))] <- paste0(100*diff(interval), '%')
    plot_data_summary <- plot_data_summary[,colnames(plot_data)]
    
    # Filter out unneeded samples
    plot_data <- subset(plot_data, s %in% sample_inds)
    
    # Combine summaries with samples
    plot_data <- rbind(plot_data, plot_data_summary)
    plot_data$Type <- factor(plot_data$Type, levels = c('Samples', 'Mean', paste0(100*diff(interval), '%')))

    plot <- ggplot(plot_data, aes(.data$Input, .data$Output, linetype = .data$Type, linewidth = as.factor(.data$s), col = .data$Type)) + 
      geom_line() +
      scale_linetype_manual(values = c(1, 1, 2)) +
      scale_linewidth_manual(values = c(rep(0.5,length(sample_inds)), 1, 1, 1), guide = 'none') +
      scale_colour_manual(values = c('grey', 'darkred', 'darkred')) +
      theme_bw() +
      labs(x = input_names[1], y = output_name, col = '', linetype = '')
    
    # If plotting multiple runs, facet by run
    if (length(inds) > 1){
      plot <- plot + facet_wrap(vars(.data$Run))
    }
    
    # If provided with true output to overlay
    if (!(is.null(Truth))){
      truth_data <- data.frame(Input = input_values, 
                               Output = c(Truth), 
                               s = max(plot_data$s)-2,
                               Run = rep(inds, each = length(output_inds)),
                               Type = 'Truth')

      plot <- plot +
        geom_line(data = truth_data) +
        scale_linetype_manual(values = c(1, 1, 2, 1)) +
        scale_colour_manual(values = c('grey', 'darkred', 'darkred', cols_good))
    }
    
    # If Obs provided, add
    if (!(is.null(Obs))){
      obs_data <- data.frame(Input = input_values[,1], 
                             Output = c(Obs), 
                             s = max(plot_data$s)-2,
                             Run = rep(inds, each = length(output_inds)),
                             Type = 'Observation')
      plot <- plot + geom_line(data = obs_data, aes(x = .data$Input, y = .data$Output), col = cols_validate[3], size = 1.5)
    }
  }
  
  
  else if (length(input_names) == 2 & !multi_var){

  }
  
  else if (length(input_names) == 2 & multi_var){
    plot_data <- data.frame(Input = input_values[,1],
                            Variable = input_values[,2],
                            Output = c(Samples), 
                            s = rep(1:ns, each = length(output_inds)),
                            Run = rep(inds, each = length(output_inds)*ns),
                            Type = 'Samples')
    
    plot_data_summary <- rbind(data.frame(stats::aggregate(Output ~ Input + Variable + Run, plot_data, mean), Type = 'Mean'),
                               data.frame(stats::aggregate(Output ~ Input + Variable + Run, plot_data, stats::quantile, probs = interval[1]), Type = 'Lower'),
                               data.frame(stats::aggregate(Output ~ Input + Variable + Run, plot_data, stats::quantile, probs = interval[2]), Type = 'Upper'))
    
    plot_data_summary$s <- max(plot_data$s)+1
    plot_data_summary$s[which(plot_data_summary$Type == 'Lower')] <- max(plot_data$s)+2 # to allow to set different linetype
    plot_data_summary$s[which(plot_data_summary$Type == 'Upper')] <- max(plot_data$s)+3 # to allow to set different linetype
    plot_data_summary$Type[which(plot_data_summary$Type %in% c('Lower','Upper'))] <- paste0(100*diff(interval), '%')
    plot_data_summary <- plot_data_summary[,colnames(plot_data)]
    
    # Filter out unneeded samples
    plot_data <- subset(plot_data, s %in% sample_inds)
    
    # Combine summaries with samples
    plot_data <- rbind(plot_data, plot_data_summary)
    plot_data$Type <- factor(plot_data$Type, levels = c('Samples', 'Mean', paste0(100*diff(interval), '%')))
    
    plot <- ggplot(plot_data, aes(.data$Input, .data$Output, linetype = .data$Type, linewidth = as.factor(.data$s), col = .data$Type)) + 
      geom_line() +
      scale_linetype_manual(values = c(1, 1, 2)) +
      scale_linewidth_manual(values = c(rep(0.5,length(sample_inds)), 1, 1, 1), guide = 'none') +
      scale_colour_manual(values = c('grey', 'darkred', 'darkred')) +
      theme_bw() +
      labs(x = input_names[1], y = output_name, col = '', linetype = '')
    
    # If plotting multiple runs, facet by run
    if (length(inds) > 1){
      plot <- plot + facet_grid(.data$Variable ~ .data$Run)
    }
    
    else {
      plot <- plot + facet_wrap(vars(.data$Variable))
    }
    
    # If provided with true output to overlay
    if (!(is.null(Truth))){
      truth_data <- data.frame(Input = input_values[,1], 
                               Variable = input_values[,2], 
                               Output = c(Truth), 
                               s = max(plot_data$s)-2,
                               Run = rep(inds, each = length(output_inds)),
                               Type = 'Truth')

      plot <- plot +
        geom_line(data = truth_data) +
        scale_linetype_manual(values = c(1, 1, 2, 1)) +
        scale_colour_manual(values = c('grey', 'darkred', 'darkred', cols_good))
    }
    
    if (!(is.null(Obs))){
      obs_data <- data.frame(Input = input_values[,1], 
                             Variable = input_values[,2], 
                             Output = c(Obs), 
                             s = max(plot_data$s)-2,
                             Run = rep(inds, each = length(output_inds)),
                             Type = 'Observation')
      
      plot <- plot + geom_line(data = obs_data, aes(x = .data$Input, y = .data$Output), linetype = 1, col = cols_validate[3], size = 1.5)
    }
    
    
  }
  
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
#' @returns description
#' 
#' @export
PlotPair <- function(coeffs, x = 'C1', y = 'C2', col = NULL, obs = NULL){
  
  coeffs <- as.data.frame(coeffs)
  plot_data <- data.frame(x = coeffs[,x],
                          y = coeffs[,y])
  
  if (!(is.null(col))){
    plot_data$col <- coeffs[,col]
  }
  
  else {
    plot_data$col <- NA
  }
  
  if (!(is.null(obs))){
    obs <- as.data.frame(obs)
    plot_data_obs <- data.frame(x = obs[,x],
                                y = obs[,y])
  }
  
  plot <- ggplot(plot_data, aes(.data$x, .data$y, col = .data$col)) +
    geom_point(size = 2) +
    labs(x = x, y = y, col = col) +
    theme_bw()
  
  if (is.null(col)){
    plot <- plot + scale_colour_discrete(guide = 'none')
  }
  
  else {
    plot <- plot + viridis::scale_colour_viridis()
  }
  
  if (!(is.null(obs))){
    plot <- plot + geom_point(data = plot_data_obs, col = cols_good, size = 4, shape = 17)
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
#' @returns description
#' 
#' @export
PlotReconError <- function(DataBasis, obs, qmax = NULL, ...){
  
  n <- ncol(DataBasis$tBasis)
  
  if (is.null(qmax)){
    qmax <- n
  }
  
  # Calculate R_W for each truncation
  RW <- errors(DataBasis$tBasis[,1:qmax], obs, ...)
  
  # Theoretical bound
  bound <- stats::qchisq(0.995, length(obs))/length(obs)
  
  plot <- ggplot(data.frame(q = 1:qmax, y = RW), aes(x = .data$q, y = .data$y)) +
    geom_line(col = 'red') +
    geom_point(col = 'red', size = 0.75) +
    ylim(0,RW[1]) +
    labs(x = 'Vector', y = 'Reconstruction Error') +
    theme_bw() +
    geom_hline(yintercept = bound, linetype = 'dashed')
  
  return(plot)
}


#' Plots proportion of variance cumulatively, or individually, explained by each basis vector
#'
#' @param DataBasis 
#' @param type 
#' @param ... 
#' 
#' @returns description 
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







#' Plotting emulator predictions
#' 
#' Plots 2D surface of emulator mean and variance across input space
#'
#' @param predictions 
#' @param design If `predictions` is the output of [Predict()], or similarly has only components `$Expectation` and `$Variance`, the corresponding set of inputs is required. Defaults to `NULL`, in which case it searches for a design in `predictions`.
#' @param inputs Names of inputs to plot on the x and y axis respectively. Defaults to `NULL`, in which case the first two columns of the design are used.
#' @param output_ind Which column of predictions to use. Defaults to `NULL`, in which case the first output is plotted.
#' @param plots Which plots to produce. Options are `mean`, `var` and `both` (the default, which plots both the mean and variance).
#'
#' @returns Plot of mean and/or variance
#' @export
#'
#' @examples
PlotPredictions <- function(predictions, design = NULL, inputs = NULL, output_ind = NULL, plots = 'both'){
  
  # Check whether predictions is output of HistoryMatch
  if (!(is.null(predictions$Preds))){
    design <- predictions$Design
    predictions <- predictions$Preds
  }
  
  # Otherwise assumes that predictions are output of Predict, and have $Expectation and $Variance
  if (is.null(predictions$Expectation)){
    stop('Expected predictions to have component $Expectation or $Preds$Expectation')
  }
  
  if (is.null(predictions$Variance)){
    stop('Expected predictions to have component $Variance or $Preds$Variance')
  }
  
  if (is.null(inputs)){
    inputs <- colnames(design)[1:2]
  }
  
  else if (!(length(inputs) == 2)){
    stop('inputs should contain exactly 2 named inputs')
  }
  
  if (!(all(inputs %in% colnames(design)))){
    stop('inputs should be contained in design or predictions$Design')
  }
  
  if (is.null(output_ind)){
    output_ind <- 1
  }
  
  if (output_ind > ncol(predictions$Expectation)){
    stop('Chosen output index too large. Please select a valid output')
  }
  
  plot_data <- data.frame(design[,inputs], 
                          Mean = predictions$Expectation[,output_ind], 
                          Variance = predictions$Variance[,output_ind])
  colnames(plot_data)[1:2] <- c('x', 'y')
  
  if (plots %in% c('both', 'mean')){
    plot1 <- ggplot(plot_data, aes(x, y, col = Mean)) +
      geom_point() +
      scale_colour_viridis_c() +
      theme_bw() +
      labs(x = inputs[1], y = inputs[2])
  }
  
  if (plots %in% c('both', 'var')){
    plot2 <- ggplot(plot_data, aes(x, y, col = Variance)) +
      geom_point() +
      scale_colour_viridis_c(option = 'A') +
      theme_bw() +
      labs(x = inputs[1], y = inputs[2])
  }
  
  if (plots == 'both'){
    plot <- cowplot::plot_grid(plot1, plot2)
  }
  
  else if (plots == 'mean'){
    plot <- plot1
  }
  
  else if (plots == 'var'){
    plot <- plot2
  }

  return(plot)
}


