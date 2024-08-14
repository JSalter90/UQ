# Load in predictions generated from 3_RunningExperiments, compare with different metrics
setwd("~/Dropbox/Exeter/Projects/UQ4Covid/CountBasis_paper")
source('code/1_ProcessData.R')
setwd("~/Dropbox/Exeter/Projects/UQ4Covid/CountBasis_paper/data/samples_final_LAD")

# Load in each in turn, store predictions
set.seed(3929)
experiments <- data.frame(n = rep(c(100,150,200), each = 50),
                          seed = sample(1:10^6, 150))

# Define columns for each total metric
experiments$RMSE_Total <- experiments$In95_Total <- NA
experiments$RMSE_Total_Het <- experiments$In95_Total_Het <- NA

for (i in 1:nrow(experiments)){
  All_em <- readRDS(paste0('Total/All_em_', experiments$seed[i], '.rds'))
  experiments$RMSE_Total[i] <- Metrics::rmse(All_em$Truth, All_em$Mean)
  experiments$In95_Total[i] <- sum(All_em$In95) / nrow(All_em)
  
  All_em_het <- readRDS(paste0('Total/All_em_het_', experiments$seed[i], '.rds'))
  experiments$RMSE_Total_Het[i] <- Metrics::rmse(All_em_het$Truth, All_em_het$Mean)
  experiments$In95_Total_Het[i] <- sum(All_em_het$In95) / nrow(All_em_het)
}

summary(experiments$RMSE_Total)
summary(experiments$RMSE_Total_Het)

summary(experiments$In95_Total)
summary(experiments$In95_Total_Het)

# Each region metric
results_all <- experiments
results_all$Region <- 'All'

for (r in 1:9){
  results_region <- results_all[1:nrow(experiments),]
  results_region$Region <- paste0('E1200000', r)
  results_region$In95_Total <- results_region$RMSE_Total <- results_region$In95_Total_Het <- results_region$RMSE_Total_Het <- NA
  for (i in 1:nrow(experiments)){
    All_em <- readRDS(paste0('E1200000', r, '/E1200000', r, '_em_', experiments$seed[i], '.rds'))
    results_region$RMSE_Total[i] <- Metrics::rmse(All_em$Truth, All_em$Mean)
    results_region$In95_Total[i] <- sum(All_em$In95) / nrow(All_em)
    
    All_em_het <- readRDS(paste0('E1200000', r, '/E1200000', r, '_em_het_', experiments$seed[i], '.rds'))
    results_region$RMSE_Total_Het[i] <- Metrics::rmse(All_em_het$Truth, All_em_het$Mean)
    results_region$In95_Total_Het[i] <- sum(All_em_het$In95) / nrow(All_em_het)
  }
  
  results_all <- rbind(results_all, results_region)
}

results_all$Region <- as.factor(results_all$Region)
summary(results_all$Region)
summary(results_all$RMSE_Total)
summary(results_all$RMSE_Total_Het)
summary(results_all$In95_Total)
summary(results_all$In95_Total_Het)

# Add equivalent metrics for each basis approach
results_all$Basis <- 'None'

# The outputs here are samples across all LADs/wards
# To compare, need to aggregate to appropriate region, and link id to true output
# First do overall total
totals <- aggregate(deaths ~ output + replicate, subset(output_LAD, week == 12 & replicate == 1), sum)
totals$LogDeaths <- log(totals$deaths)

results_tmp <- results_all[1:nrow(experiments),]
results_tmp$Region <- 'All'
results_tmp$Basis <- 'PLNPCA'
results_tmp$In95_Total <- results_tmp$RMSE_Total <- results_tmp$In95_Total_Het <- results_tmp$RMSE_Total_Het <- NA

for (i in 1:nrow(experiments)){
  # Check that simulation exists
  test <- list.files('PLNPCA/', pattern = paste0(experiments$seed[i]))
  if (length(test) > 0){
    All_em <- readRDS(paste0('PLNPCA/basis_PLNPCA_em_', experiments$seed[i], '.rds'))
    tmp <- AggregateSamples(All_em$samples)
    tmp2 <- data.frame(Mean = log(tmp$mean),
                       Lower = log(tmp$lower),
                       Upper = log(tmp$upper),
                       Output = All_em$Run)
    tmp2$Truth <- totals$LogDeaths[match(All_em$Run, totals$output)]
    tmp2$In95 <- tmp2$Truth >= tmp2$Lower & tmp2$Truth <= tmp2$Upper
    results_tmp$RMSE_Total[i] <- Metrics::rmse(tmp2$Truth, tmp2$Mean)
    results_tmp$In95_Total[i] <- sum(tmp2$In95) / nrow(tmp2)
    
    All_em <- readRDS(paste0('PLNPCA/basis_PLNPCA_em_het_', experiments$seed[i], '.rds'))
    tmp <- AggregateSamples(All_em$samples)
    tmp2 <- data.frame(Mean = log(tmp$mean),
                       Lower = log(tmp$lower),
                       Upper = log(tmp$upper),
                       Output = All_em$Run)
    tmp2$Truth <- totals$LogDeaths[match(All_em$Run, totals$output)]
    tmp2$In95 <- tmp2$Truth >= tmp2$Lower & tmp2$Truth <= tmp2$Upper
    results_tmp$RMSE_Total_Het[i] <- Metrics::rmse(tmp2$Truth, tmp2$Mean)
    results_tmp$In95_Total_Het[i] <- sum(tmp2$In95) / nrow(tmp2)
  }
}

results_all <- rbind(results_all, results_tmp)

# Repeat by region
# Need to be careful that compare properly when did log(x+1) for the individual emulator
totals <- aggregate(deaths ~ output + replicate + region, subset(output_LAD, week == 12 & replicate == 1), sum)
totals$LogDeaths <- log(totals$deaths)
totals$LogDeaths1 <- log(totals$deaths + 1)

# Combine indices for each region into a list
ind_region <- list(LAD_NE, LAD_NW, LAD_YO, LAD_EM, LAD_WM, LAD_EE, LAD_LO, LAD_SE, LAD_SW)

# List all LAD indices
location_ids <- unique(output_LAD$LAD19CD)

for (r in 1:9){
  results_region <- results_all[1:nrow(experiments),]
  results_region$Region <- paste0('E1200000', r)
  results_region$Basis <- 'PLNPCA'
  results_region$In95_Total <- results_region$RMSE_Total <- results_region$In95_Total_Het <- results_region$RMSE_Total_Het <- NA
  
  totals_region <- subset(totals, region == paste0('E1200000', r))
  
  for (i in 1:nrow(experiments)){
    # Check that simulation exists
    test <- list.files('PLNPCA/', pattern = paste0(experiments$seed[i]))
    if (length(test) > 0){
      All_em <- readRDS(paste0('PLNPCA/basis_PLNPCA_em_', experiments$seed[i], '.rds'))
      tmp <- AggregateSamples(All_em$samples, locs = which(location_ids %in% ind_region[[r]]))
      if (any(totals_region$deaths == 0)){
        tmp2 <- data.frame(Mean = log(tmp$mean + 1),
                           Lower = log(tmp$lower + 1),
                           Upper = log(tmp$upper + 1),
                           Output = All_em$Run)
        tmp2$Truth <- totals_region$LogDeaths1[match(All_em$Run, totals_region$output)]
      }
      else {
        tmp2 <- data.frame(Mean = log(tmp$mean),
                           Lower = log(tmp$lower),
                           Upper = log(tmp$upper),
                           Output = All_em$Run)
        tmp2$Truth <- totals_region$LogDeaths[match(All_em$Run, totals_region$output)]
      }
      tmp2$In95 <- tmp2$Truth >= tmp2$Lower & tmp2$Truth <= tmp2$Upper
      results_region$RMSE_Total[i] <- Metrics::rmse(tmp2$Truth, tmp2$Mean)
      results_region$In95_Total[i] <- sum(tmp2$In95) / nrow(tmp2)
      
      All_em <- readRDS(paste0('PLNPCA/basis_PLNPCA_em_het_', experiments$seed[i], '.rds'))
      tmp <- AggregateSamples(All_em$samples, locs = which(location_ids %in% ind_region[[r]]))
      if (any(totals_region$deaths == 0)){
        tmp2 <- data.frame(Mean = log(tmp$mean + 1),
                           Lower = log(tmp$lower + 1),
                           Upper = log(tmp$upper + 1),
                           Output = All_em$Run)
        tmp2$Truth <- totals_region$LogDeaths1[match(All_em$Run, totals_region$output)]
      }
      else {
        tmp2 <- data.frame(Mean = log(tmp$mean),
                           Lower = log(tmp$lower),
                           Upper = log(tmp$upper),
                           Output = All_em$Run)
        tmp2$Truth <- totals_region$LogDeaths[match(All_em$Run, totals_region$output)]
      }
      tmp2$In95 <- tmp2$Truth >= tmp2$Lower & tmp2$Truth <= tmp2$Upper
      results_region$RMSE_Total_Het[i] <- Metrics::rmse(tmp2$Truth, tmp2$Mean)
      results_region$In95_Total_Het[i] <- sum(tmp2$In95) / nrow(tmp2)
    }
  }
  results_all <- rbind(results_all, results_region)
}

results_all$Basis <- as.factor(results_all$Basis)

saveRDS(results_all, file = 'results_all.rds')



# Repeating for SVD
results_all <- readRDS('results_all.rds')

totals <- aggregate(deaths ~ output + replicate, subset(output_LAD, week == 12 & replicate == 1), sum)
totals$LogDeaths <- log(totals$deaths)

results_tmp <- results_all[1:nrow(experiments),]
results_tmp$Region <- 'All'
results_tmp$Basis <- 'SVD'
results_tmp$In95_Total <- results_tmp$RMSE_Total <- results_tmp$In95_Total_Het <- results_tmp$RMSE_Total_Het <- NA

for (i in 1:nrow(experiments)){
  # Check that simulation exists
  test <- list.files('SVD/', pattern = paste0(experiments$seed[i]))
  if (length(test) > 0){
    All_em <- readRDS(paste0('SVD/basis_SVD_em_', experiments$seed[i], '.rds'))
    tmp <- AggregateSamples(All_em$samples)
    tmp2 <- data.frame(Mean = log(tmp$mean),
                       Lower = log(tmp$lower),
                       Upper = log(tmp$upper),
                       Output = All_em$Run)
    tmp2$Truth <- totals$LogDeaths[match(All_em$Run, totals$output)]
    tmp2$In95 <- tmp2$Truth >= tmp2$Lower & tmp2$Truth <= tmp2$Upper
    results_tmp$RMSE_Total[i] <- Metrics::rmse(tmp2$Truth, tmp2$Mean)
    results_tmp$In95_Total[i] <- sum(tmp2$In95) / nrow(tmp2)
    
    All_em <- readRDS(paste0('SVD/basis_SVD_em_het_', experiments$seed[i], '.rds'))
    tmp <- AggregateSamples(All_em$samples)
    tmp2 <- data.frame(Mean = log(tmp$mean),
                       Lower = log(tmp$lower),
                       Upper = log(tmp$upper),
                       Output = All_em$Run)
    tmp2$Truth <- totals$LogDeaths[match(All_em$Run, totals$output)]
    tmp2$In95 <- tmp2$Truth >= tmp2$Lower & tmp2$Truth <= tmp2$Upper
    results_tmp$RMSE_Total_Het[i] <- Metrics::rmse(tmp2$Truth, tmp2$Mean)
    results_tmp$In95_Total_Het[i] <- sum(tmp2$In95) / nrow(tmp2)
  }
}

results_all <- rbind(results_all, results_tmp)

# Repeat by region
totals <- aggregate(deaths ~ output + replicate + region, subset(output_LAD, week == 12 & replicate == 1), sum)
totals$LogDeaths <- log(totals$deaths)
totals$LogDeaths1 <- log(totals$deaths + 1)

# Combine indices for each region into a list
ind_region <- list(LAD_NE, LAD_NW, LAD_YO, LAD_EM, LAD_WM, LAD_EE, LAD_LO, LAD_SE, LAD_SW)

# List all LAD indices
location_ids <- unique(output_LAD$LAD19CD)

for (r in 1:9){
  results_region <- results_all[1:nrow(experiments),]
  results_region$Region <- paste0('E1200000', r)
  results_region$Basis <- 'SVD'
  results_region$In95_Total <- results_region$RMSE_Total <- results_region$In95_Total_Het <- results_region$RMSE_Total_Het <- NA
  
  totals_region <- subset(totals, region == paste0('E1200000', r))
  
  for (i in 1:nrow(experiments)){
    # Check that simulation exists
    test <- list.files('SVD/', pattern = paste0(experiments$seed[i]))
    if (length(test) > 0){
      All_em <- readRDS(paste0('SVD/basis_SVD_em_', experiments$seed[i], '.rds'))
      tmp <- AggregateSamples(All_em$samples, locs = which(location_ids %in% ind_region[[r]]))
      if (any(totals_region$deaths == 0)){
        tmp2 <- data.frame(Mean = log(tmp$mean + 1),
                           Lower = log(tmp$lower + 1),
                           Upper = log(tmp$upper + 1),
                           Output = All_em$Run)
        tmp2$Truth <- totals_region$LogDeaths1[match(All_em$Run, totals_region$output)]
      }
      else {
        tmp2 <- data.frame(Mean = log(tmp$mean),
                           Lower = log(tmp$lower),
                           Upper = log(tmp$upper),
                           Output = All_em$Run)
        tmp2$Truth <- totals_region$LogDeaths[match(All_em$Run, totals_region$output)]
      }
      tmp2$In95 <- tmp2$Truth >= tmp2$Lower & tmp2$Truth <= tmp2$Upper
      results_region$RMSE_Total[i] <- Metrics::rmse(tmp2$Truth, tmp2$Mean)
      results_region$In95_Total[i] <- sum(tmp2$In95) / nrow(tmp2)
      
      All_em <- readRDS(paste0('SVD/basis_SVD_em_het_', experiments$seed[i], '.rds'))
      tmp <- AggregateSamples(All_em$samples, locs = which(location_ids %in% ind_region[[r]]))
      if (any(totals_region$deaths == 0)){
        tmp2 <- data.frame(Mean = log(tmp$mean + 1),
                           Lower = log(tmp$lower + 1),
                           Upper = log(tmp$upper + 1),
                           Output = All_em$Run)
        tmp2$Truth <- totals_region$LogDeaths1[match(All_em$Run, totals_region$output)]
      }
      else {
        tmp2 <- data.frame(Mean = log(tmp$mean),
                           Lower = log(tmp$lower),
                           Upper = log(tmp$upper),
                           Output = All_em$Run)
        tmp2$Truth <- totals_region$LogDeaths[match(All_em$Run, totals_region$output)]
      }
      tmp2$In95 <- tmp2$Truth >= tmp2$Lower & tmp2$Truth <= tmp2$Upper
      results_region$RMSE_Total_Het[i] <- Metrics::rmse(tmp2$Truth, tmp2$Mean)
      results_region$In95_Total_Het[i] <- sum(tmp2$In95) / nrow(tmp2)
    }
  }
  results_all <- rbind(results_all, results_region)
}

results_all$Basis <- as.factor(results_all$Basis)

saveRDS(results_all, file = 'results_all.rds')






#### Repeating for ward-level samples ####
# Larger files, stored separately
# Create results file from scratch, then combine with results_all
setwd("/Volumes/Extreme_SSD/data/CountBasis/samples_final_ward")
set.seed(3929)
experiments <- data.frame(n = rep(c(100,150,200), each = 50),
                          seed = sample(1:10^6, 150))
results_ward <- experiments
results_ward$RMSE_Total <- results_ward$In95_Total <- NA
results_ward$RMSE_Total_Het <- results_ward$In95_Total_Het <- NA
results_ward$Region <- 'All'
results_ward$Basis <- 'PLNPCA_ward'

totals <- aggregate(deaths ~ output + replicate, subset(output_ward, replicate == 1), sum)
totals$LogDeaths <- log(totals$deaths)

for (i in 1:nrow(experiments)){
  # Check that simulation exists
  test <- list.files('PLNPCA/', pattern = paste0(experiments$seed[i]))
  if (length(test) > 0){
    All_em <- readRDS(paste0('PLNPCA/basis_PLNPCA_em_', experiments$seed[i], '.rds'))
    tmp <- AggregateSamples(All_em$samples)
    tmp2 <- data.frame(Mean = log(tmp$mean),
                       Lower = log(tmp$lower),
                       Upper = log(tmp$upper),
                       Output = All_em$Run)
    tmp2$Truth <- totals$LogDeaths[match(All_em$Run, totals$output)]
    tmp2$In95 <- tmp2$Truth >= tmp2$Lower & tmp2$Truth <= tmp2$Upper
    results_ward$RMSE_Total[i] <- Metrics::rmse(tmp2$Truth, tmp2$Mean)
    results_ward$In95_Total[i] <- sum(tmp2$In95) / nrow(tmp2)
    
    All_em <- readRDS(paste0('PLNPCA/basis_PLNPCA_em_het_', experiments$seed[i], '.rds'))
    tmp <- AggregateSamples(All_em$samples)
    tmp2 <- data.frame(Mean = log(tmp$mean),
                       Lower = log(tmp$lower),
                       Upper = log(tmp$upper),
                       Output = All_em$Run)
    tmp2$Truth <- totals$LogDeaths[match(All_em$Run, totals$output)]
    tmp2$In95 <- tmp2$Truth >= tmp2$Lower & tmp2$Truth <= tmp2$Upper
    results_ward$RMSE_Total_Het[i] <- Metrics::rmse(tmp2$Truth, tmp2$Mean)
    results_ward$In95_Total_Het[i] <- sum(tmp2$In95) / nrow(tmp2)
  }
}

setwd("~/Dropbox/Exeter/Projects/UQ4Covid/CountBasis_paper")
saveRDS(results_ward, file = 'data/results_ward.rds')

# Repeat by region
totals <- aggregate(deaths ~ output + replicate + region, subset(output_ward, replicate == 1), sum)
totals$LogDeaths <- log(totals$deaths)
totals$LogDeaths1 <- log(totals$deaths + 1)

# Combine indices for each region into a list
ind_region <- list(ward_NE, ward_NW, ward_YO, ward_EM, ward_WM, ward_EE, ward_LO, ward_SE, ward_SW)

# List all LAD indices
location_ids <- unique(output_ward$ward)
setwd("/Volumes/Extreme_SSD/data/CountBasis/samples_final_ward")
for (r in 1:9){
  results_region <- results_ward[1:nrow(experiments),]
  results_region$Region <- paste0('E1200000', r)
  results_region$Basis <- 'PLNPCA_ward'
  results_region$In95_Total <- results_region$RMSE_Total <- results_region$In95_Total_Het <- results_region$RMSE_Total_Het <- NA
  
  totals_region <- subset(totals, region == paste0('E1200000', r))
  
  for (i in 1:nrow(experiments)){
    # Check that simulation exists
    test <- list.files('PLNPCA/', pattern = paste0(experiments$seed[i]))
    if (length(test) > 0){
      All_em <- readRDS(paste0('PLNPCA/basis_PLNPCA_em_', experiments$seed[i], '.rds'))
      
      if (dim(All_em$samples)[1] == 8071){
        location_ids <- unique(output_ward$ward)
      }
      if (dim(All_em$samples)[1] == 8063){
        location_ids <- unique(subset(output_ward, !(ward %in% c("E05011090","E05011092","E05011094","E05009493",
                                                                 "E05011788","E05011924","E05010147","E05012767")))$ward)
      }
      if (dim(All_em$samples)[1] == 8060){
        location_ids <- unique(subset(output_ward, !(ward %in% c("E05011090","E05011092","E05011094","E05009493",
                                                                 "E05011788","E05011924","E05010147","E05012767",
                                                                 "E05008329","E05012716","E05011898")))$ward)
      }
      if (experiments$seed[i] == 45138){
        location_ids <- unique(subset(output_ward, !(ward %in% c("E05011090","E05011092","E05011094","E05009493",
                                                                 "E05011788","E05011924","E05010147","E05012767",
                                                                 "E05008329","E05012716","E05011898",
                                                                 "E05011091","E05011644")))$ward)
      }
      if (experiments$seed[i] == 232487){
        location_ids <- unique(subset(output_ward, !(ward %in% c("E05011090","E05011092","E05011094","E05009493",
                                                                 "E05011788","E05011924","E05010147","E05012767",
                                                                 "E05008329","E05012716","E05011898",
                                                                 "E05008330","E05011201","E05012762")))$ward)
      }
      if (experiments$seed[i] == 945202){
        location_ids <- unique(subset(output_ward, !(ward %in% c("E05011090","E05011092","E05011094","E05009493",
                                                                 "E05011788","E05011924","E05010147","E05012767",
                                                                 "E05008329","E05012716","E05011898",
                                                                 "E05010674","E05010061")))$ward)
      }
      tmp <- AggregateSamples(All_em$samples, locs = which(location_ids %in% ind_region[[r]]))
      if (any(totals_region$deaths == 0)){
        tmp2 <- data.frame(Mean = log(tmp$mean + 1),
                           Lower = log(tmp$lower + 1),
                           Upper = log(tmp$upper + 1),
                           Output = All_em$Run)
        tmp2$Truth <- totals_region$LogDeaths1[match(All_em$Run, totals_region$output)]
      }
      else {
        tmp2 <- data.frame(Mean = log(tmp$mean),
                           Lower = log(tmp$lower),
                           Upper = log(tmp$upper),
                           Output = All_em$Run)
        tmp2$Truth <- totals_region$LogDeaths[match(All_em$Run, totals_region$output)]
      }
      tmp2$In95 <- tmp2$Truth >= tmp2$Lower & tmp2$Truth <= tmp2$Upper
      results_region$RMSE_Total[i] <- Metrics::rmse(tmp2$Truth, tmp2$Mean)
      results_region$In95_Total[i] <- sum(tmp2$In95) / nrow(tmp2)
      
      All_em <- readRDS(paste0('PLNPCA/basis_PLNPCA_em_het_', experiments$seed[i], '.rds'))
      tmp <- AggregateSamples(All_em$samples, locs = which(location_ids %in% ind_region[[r]]))
      if (any(totals_region$deaths == 0)){
        tmp2 <- data.frame(Mean = log(tmp$mean + 1),
                           Lower = log(tmp$lower + 1),
                           Upper = log(tmp$upper + 1),
                           Output = All_em$Run)
        tmp2$Truth <- totals_region$LogDeaths1[match(All_em$Run, totals_region$output)]
      }
      else {
        tmp2 <- data.frame(Mean = log(tmp$mean),
                           Lower = log(tmp$lower),
                           Upper = log(tmp$upper),
                           Output = All_em$Run)
        tmp2$Truth <- totals_region$LogDeaths[match(All_em$Run, totals_region$output)]
      }
      tmp2$In95 <- tmp2$Truth >= tmp2$Lower & tmp2$Truth <= tmp2$Upper
      results_region$RMSE_Total_Het[i] <- Metrics::rmse(tmp2$Truth, tmp2$Mean)
      results_region$In95_Total_Het[i] <- sum(tmp2$In95) / nrow(tmp2)
    }
  }
  results_ward <- rbind(results_ward, results_region)
}

#results_all$Basis <- as.factor(results_all$Basis)
setwd("~/Dropbox/Exeter/Projects/UQ4Covid/CountBasis_paper")
saveRDS(results_ward, file = 'data/results_ward.rds')





#### Summarising results for tables ####
library(dplyr)
setwd("~/Dropbox/Exeter/Projects/UQ4Covid/CountBasis_paper")
results_all <- readRDS('data/samples_final_LAD/results_all.rds')
results_ward <- readRDS('data/results_ward.rds')
results_ward$Region <- as.factor(results_ward$Region)

# Select only the HetGP results, except for aggregations (Basis = 'None')
tmp <- subset(results_all, Basis == 'None')
tmp2 <- subset(results_all[,-c(3:4)], Basis == 'PLNPCA')
tmp3 <- subset(results_all[,-c(3:4)], Basis == 'SVD')
tmp$Basis <- tmp2$Basis <- tmp3$Basis <- NULL
colnames(tmp)[3:6] <- c('In95_Single_GP', 'RMSE_Single_GP', 'In95_Single', 'RMSE_Single')
colnames(tmp2)[3:4] <- c('In95_Basis', 'RMSE_Basis')
colnames(tmp3)[3:4] <- c('In95_SVD', 'RMSE_SVD')
tmp <- tmp %>% left_join(tmp2, by = c('n', 'seed', 'Region'))
tmp <- tmp %>% left_join(tmp3, by = c('n', 'seed', 'Region'))

# Add in ward info
tmp4 <- subset(results_ward[,-c(3:4)], Basis == 'PLNPCA_ward')
tmp4$Basis <- NULL
colnames(tmp4)[3:4] <- c('In95_Basis_Ward', 'RMSE_Basis_Ward')
tmp <- tmp %>% left_join(tmp4, by = c('n', 'seed', 'Region'))

# Change region codes to names
levels(tmp$Region) <- c('All', 'North East', 'North West', 'Yorkshire',
                        'East Midlands', 'West Midlands', 'East of England',
                        'London', 'South East', 'South West')

# % of experiments that are better than the baseline (HetGP + single emulator)
tmp$BetterGP <- tmp$RMSE_Single_GP < tmp$RMSE_Single
tmp$BetterLAD <- tmp$RMSE_Basis < tmp$RMSE_Single
tmp$BetterWard <- tmp$RMSE_Basis_Ward < tmp$RMSE_Single
tmp$BetterSVD <- tmp$RMSE_SVD < tmp$RMSE_Single
aggregate(cbind(BetterGP, BetterLAD, BetterWard, BetterSVD) ~ Region, data = tmp, mean)
aggregate(cbind(BetterGP, BetterLAD, BetterWard, BetterSVD) ~ Region, data = subset(tmp, n == 200), mean)
aggregate(cbind(BetterGP, BetterLAD, BetterWard, BetterSVD) ~ Region, data = subset(tmp, n == 150), mean)
aggregate(cbind(BetterGP, BetterLAD, BetterWard, BetterSVD) ~ Region, data = subset(tmp, n == 100), mean)

# Median % in 95% interval
aggregate(cbind(In95_Single,In95_Single_GP,In95_SVD,In95_Basis,In95_Basis_Ward) ~ Region, data = tmp, median)
aggregate(cbind(In95_Single,In95_Single_GP,In95_SVD,In95_Basis,In95_Basis_Ward) ~ Region, data = subset(tmp, n == 200), median)
aggregate(cbind(In95_Single,In95_Single_GP,In95_SVD,In95_Basis,In95_Basis_Ward) ~ Region, data = subset(tmp, n == 150), median)
aggregate(cbind(In95_Single,In95_Single_GP,In95_SVD,In95_Basis,In95_Basis_Ward) ~ Region, data = subset(tmp, n == 100), median)

# % change in RMSE relative to baseline
# Find median by region
results_summary <- aggregate(cbind(RMSE_Single,RMSE_Single_GP,RMSE_SVD,RMSE_Basis,RMSE_Basis_Ward) ~ Region, data = subset(tmp, n == 200), median)
results_summary$ChangeGP <- results_summary$RMSE_Single_GP / results_summary$RMSE_Single
results_summary$ChangeLAD <- results_summary$RMSE_Basis / results_summary$RMSE_Single
results_summary$ChangeWard <- results_summary$RMSE_Basis_Ward / results_summary$RMSE_Single
results_summary$ChangeSVD <- results_summary$RMSE_SVD / results_summary$RMSE_Single
results_summary

results_summary <- aggregate(cbind(RMSE_Single,RMSE_Single_GP,RMSE_SVD,RMSE_Basis,RMSE_Basis_Ward) ~ Region, data = subset(tmp, n == 150), median)
results_summary$ChangeGP <- results_summary$RMSE_Single_GP / results_summary$RMSE_Single
results_summary$ChangeLAD <- results_summary$RMSE_Basis / results_summary$RMSE_Single
results_summary$ChangeWard <- results_summary$RMSE_Basis_Ward / results_summary$RMSE_Single
results_summary$ChangeSVD <- results_summary$RMSE_SVD / results_summary$RMSE_Single
results_summary

results_summary <- aggregate(cbind(RMSE_Single,RMSE_Single_GP,RMSE_SVD,RMSE_Basis,RMSE_Basis_Ward) ~ Region, data = subset(tmp, n == 100), median)
results_summary$ChangeGP <- results_summary$RMSE_Single_GP / results_summary$RMSE_Single
results_summary$ChangeLAD <- results_summary$RMSE_Basis / results_summary$RMSE_Single
results_summary$ChangeWard <- results_summary$RMSE_Basis_Ward / results_summary$RMSE_Single
results_summary$ChangeSVD <- results_summary$RMSE_SVD / results_summary$RMSE_Single
results_summary

