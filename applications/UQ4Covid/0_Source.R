# Useful functions, loading packages etc.
# General functions will likely end up in code
library(ggplot2)
library(PLNmodels)

# Process into data matrix
ProcessData <- function(data, by_id = 'LAD19CD', output = 'deaths'){
  # Info on locations
  location_ids <- unique(data[,by_id])
  ell <- length(location_ids)
  
  # Info on runs
  run_ids <- unique(data[,1:2]) # need combination of run id and duplicate number
  n <- nrow(run_ids)
  
  # Is the provided dataset complete?
  stopifnot(n == nrow(data) / ell)
  
  data_matrix <- matrix(0, ell, n)
  
  for (i in 1:n){
    tmp_data <- subset(data, output == run_ids$output[i] & replicate == run_ids$replicate[i])[,c(by_id,output)]
    
    # Map to same ordering as in location_ids
    tmp_data <- left_join(as.data.frame(location_ids), 
                          tmp_data, 
                          by = c('location_ids' = paste0(by_id)))
    
    # Add to matrix
    data_matrix[,i] <- tmp_data[,output]
  }
  
  return(list(data = data_matrix,
              location = location_ids,
              run = run_ids))
}

# Example
output_LAD <- readRDS('data/output_LAD.rds')
ens_data <- ProcessData(subset(output_LAD, week == 12 & replicate == 1))

# Some checks that worked
tmp1 <- apply(ens_data$data, 2, sum) # total deaths per run
tmp2 <- aggregate(deaths ~ output, subset(output_LAD, week == 12 & replicate == 1), sum)$deaths
summary(tmp1 - tmp2)

tmp1 <- apply(ens_data$data, 1, sum) # total deaths per LAD across 250 runs
tmp2 <- aggregate(deaths ~ LAD19CD, subset(output_LAD, week == 12 & replicate == 1), sum)$deaths
summary(tmp1 - tmp2)

# Hospitalisations instead
ens_data <- ProcessData(subset(output_LAD, week == 12 & replicate == 1), output = 'cumH')
tmp1 <- apply(ens_data$data, 2, sum)
tmp2 <- aggregate(cumH ~ output, subset(output_LAD, week == 12 & replicate == 1), sum)$cumH
summary(tmp1 - tmp2)

tmp1 <- apply(ens_data$data, 1, sum)
tmp2 <- aggregate(cumH ~ LAD19CD, subset(output_LAD, week == 12 & replicate == 1), sum)$cumH
summary(tmp1 - tmp2)


# Finding basis
?PLNPCA
?prepare_data

CountBasis <- function(data){
  
  
  
}



