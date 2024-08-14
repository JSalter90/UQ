setwd("~/Dropbox/Exeter/Projects/UQ4Covid/CountBasis_paper")
source('code/1_ProcessData.R') # also loads in 0_CountBasis.R
setwd('~/Dropbox/UQ/')
source('code/PlotFunctions.R')
source('code/Gasp.R')
setwd("~/Dropbox/Exeter/Projects/UQ4Covid/CountBasis_paper")

# Inputs
dim(design) # 250 unique input sets, $repeats gives number of replicates for each
# Outputs
dim(output_LAD) # 1661100 = 339 LADs x 7 weeks x 700 simulations
dim(output_ward) # 5649700 = 8071 wards x 700 simulations

# Plot some summaries of output, e.g. total hospitalisations vs deaths in week 12
output_total <- aggregate(cbind(cumH, deaths) ~ output + replicate, data = subset(output_LAD, week == 12), sum)
ggplot(output_total, aes(log(cumH), log(deaths))) + geom_point()

# Specific example, for random split into training/validation
# As seen in Results/Worked Example section
# Just demonstrating how works here, comparison to single emulators comes later

#### LAD level ####
# Split into training/validation
n <- nrow(design)
train_prop <- 0.8
val_prop <- 1 - train_prop

set.seed(34810)
train_inds <- sample(1:n, train_prop*n)
train_labels <- design$output[train_inds]
# Just use single replicate for each input
train_data <- subset(output_LAD, week == 12 & replicate == 1 & output %in% train_labels) # 67800 = 200 runs x 339 LADs
val_data <- subset(output_LAD, week == 12 & replicate == 1 & !(output %in% train_labels)) # 16950 = 50 runs x 339 LADs

# Process data into ell x n
train_data <- ProcessData(train_data)
val_data <- ProcessData(val_data)

dim(train_data$data) # 339 x 200
dim(val_data$data) # 339 x 50

head(train_data$run) # labels corresponding to the columns of $data
head(train_data$location) # labels corresponding to rows of $data
head(val_data$location) # should be same ordering in training and validation sets

# Construct basis
# basis_deaths <- CountBasis(train_data$data, rank = 10) # rank is a design choice, use small value here for speed
# saveRDS(basis_deaths, file = 'data/basis_deaths_example.rds')
basis_deaths <- readRDS('data/basis_deaths_example.rds')
summary(basis_deaths)
dim(basis_deaths$tBasis) # 339x10

# Plot leading basis vectors
png('plots/basis_LAD.png', width = 1344, height = 960, res = 200)
PlotBasis(basis_deaths, q = 4, LonLat[match(train_data$location, LonLat$LAD2019CD),2:3])
dev.off()

# Extract coefficients, combine with design
# Need to use tags in $run$output to match up with design
tData <- GetEmDataCount(design_em[match(train_data$run$output, design_em$output),1:15], basis_deaths, q = 10)

# Plotting pairs of coefficients, inputs
PlotPair(tData, 'C1', 'C2', col = 'R0')
PlotPair(tData, 'R0', 'ns', col = 'C1')

# Now just need to emulate the latent coefficients, for the leading q basis vectors
# For comparison, try both a standard stationary homoscedastic GP, vs hetGP
EmDeaths <- BasisEmulators(tData, 6, mean_fn = 'step', training_prop = 1)
PlotActive(EmDeaths, InputNames = colnames(tData)[1:15])
LOOs <- lapply(1:length(EmDeaths), function(k) LeaveOneOut(EmDeaths[[k]]))
cowplot::plot_grid(LOOs[[1]], LOOs[[2]], LOOs[[3]], LOOs[[4]], LOOs[[5]], LOOs[[6]])

# Predicting for validation points
val_inputs <- subset(design_em, !(output %in% train_data$run$output))
Preds_val <- BasisPredGasp(val_inputs[,1:15], EmDeaths)

# We don't know the values of the latent coefficients for runs we didn't include when estimating the basis
# For validation, therefore easiest to a) predict latent coefficients b) sample latent fields c) reconstruct true fields from latent samples
# Therefore get prediction for the total deaths, or regional deaths, or any other aggregation
val_response <- CountBasisEmSamples(Preds_val, basis_deaths, ReturnAll = TRUE, BasisUncertainty = TRUE)
dim(val_response$samples) # 339 LADs x 1000 samples x 50 validation points
ValidateSum(val_response$samples, val_data$data)

# We can restrict to particular regions, e.g. NE only
ValidateSum(val_response$samples, val_data$data, locs = which(val_data$location %in% LAD_NE))

# We can also plot all the individual samples and the truth at a LAD level, for the validation runs, e.g for London region
PlotSamplesCount(val_response$samples[,1:100,], locs = which(val_data$location %in% LAD_LO), runs = 1:16, Truth = val_data$data)

# Or all 339 LADs
PlotSamplesCount(val_response$samples[,1:100,], runs = 1, Truth = val_data$data)

# Likely better for stochastic data: hetGP
EmDeathsHet <- BasisEmulatorsHet(tData, 3, training_prop = 1)
Preds_val_Het <- BasisPredHet(val_inputs[,1:15], EmDeathsHet)
val_response_Het <- CountBasisEmSamples(Preds_val_Het, basis_deaths, ReturnAll = TRUE)

# Store objects for reproducing plots later
saveRDS(val_response_Het, 'data/val_het_LAD.rds')
saveRDS(val_data, 'data/val_data_LAD.rds')
saveRDS(Preds_val_Het, 'data/preds_het_LAD.rds')

cowplot::plot_grid(ValidateSum(val_response$samples, val_data$data), 
                   ValidateSum(val_response_Het$samples, val_data$data), nrow = 1)
Metrics::rmse(apply(val_data$data,2,sum), apply(val_response$samples,3,sum))
Metrics::rmse(apply(val_data$data,2,sum), apply(val_response_Het$samples,3,sum))
Metrics::rmse(log(apply(val_data$data,2,sum)+1), log(apply(val_response$samples,3,sum))+1)
Metrics::rmse(log(apply(val_data$data,2,sum)+1), log(apply(val_response_Het$samples,3,sum))+1)

cowplot::plot_grid(ValidateSum(val_response$samples, val_data$data, locs = which(val_data$location %in% LAD_NE)), 
                   ValidateSum(val_response_Het$samples, val_data$data, locs = which(val_data$location %in% LAD_NE)), nrow = 1)

# NE samples
PlotSamplesCount(val_response_Het$samples[,1:100,], locs = which(val_data$location %in% LAD_NE), runs = 1:16, Truth = val_data$data)

# Specific locations
ValidateSum(val_response_Het$samples, val_data$data, locs = 2)




#### Including replicates ####
# Repeat the above, but now include replicates where available
# Split replicated and unreplicated design points in same proportion
set.seed(4929)
design_single <- subset(design, repeats == 1)
design_reps <- subset(design, repeats == 10)
train_inds_single <- sample(1:nrow(design_single), train_prop*nrow(design_single))
train_inds_reps <- sample(1:nrow(design_reps), train_prop*nrow(design_reps))
train_labels <- c(design_single$output[train_inds_single], design_reps$output[train_inds_reps])
train_data <- subset(output_LAD, week == 12 & output %in% train_labels)
val_data <- subset(output_LAD, week == 12 & !(output %in% train_labels))

# Process data into matrix
train_data <- ProcessData(train_data)
val_data <- ProcessData(val_data)

dim(train_data$data) # 339 x 560
dim(val_data$data) # 339 x 140
tail(train_data$run) # combinations of tag, replicate

# Basis (slower - more training data)
# basis_deaths_rep <- CountBasis(train_data$data, rank = 100)
# saveRDS(basis_deaths_rep, file = 'data/basis_deaths_rep.rds')
basis_deaths_rep <- readRDS('data/basis_deaths_rep.rds')
tData_rep <- GetEmDataCount(design_em[match(train_data$run$output, design_em$output),1:15], basis_deaths_rep, q = 10)

# To assess the stochasticity, can plot the variability in the coefficients for sets of replicates
CoeffPlots <- lapply(1:8, function(k) PlotCoeffs(basis_deaths_rep, k, run_ids = train_data$run, plot_inds = 1:12))
saveRDS(CoeffPlots, file = 'data/CoeffPlots.rds')
cowplot::plot_grid(CoeffPlots[[1]],CoeffPlots[[2]],CoeffPlots[[3]],CoeffPlots[[4]])
cowplot::plot_grid(CoeffPlots[[5]],CoeffPlots[[6]],CoeffPlots[[7]],CoeffPlots[[8]])

# hetGP should handle this dataset well
EmDeathsHet_rep <- BasisEmulatorsHet(tData_rep, 3, training_prop = 1)
val_inputs <- design_em[match(val_data$run$output, design_em$output),1:15] # including repeated inputs, to match up with output
Preds_val_Het_rep <- BasisPredHet(val_inputs[,1:15], EmDeathsHet_rep)
val_response_Het_rep <- CountBasisEmSamples(Preds_val_Het_rep, basis_deaths_rep, ReturnAll = TRUE)
ValidateSum(val_response_Het_rep$samples, val_data$data)

# Store objects for reproducing plots later
saveRDS(val_response_Het_rep, 'data/val_het_LAD_reps.rds')
saveRDS(val_data, 'data/val_data_LAD_reps.rds')

# As before, we can sample from the emulator and plot these vs different locations
# We can also plot the multiple replicates we observe vs these samples
# The uncertainty across the samples is accounting for both a) extrapolation to unobserved inputs x and b) stochasticity in output at x
k <- 1
rep_data <- ProcessData(subset(output_LAD, week == 12 & output == val_data$run$output[k]))
plot1 <- PlotReplicates(val_response_Het_rep$samples[,1:100,k], locs = which(val_data$location %in% LAD_NE), Replicates = rep_data$data)
k <- 3
rep_data <- ProcessData(subset(output_LAD, week == 12 & output == val_data$run$output[k]))
plot2 <- PlotReplicates(val_response_Het_rep$samples[,1:100,k], locs = which(val_data$location %in% LAD_NE), Replicates = rep_data$data)
k <- 5
rep_data <- ProcessData(subset(output_LAD, week == 12 & output == val_data$run$output[k]))
plot3 <- PlotReplicates(val_response_Het_rep$samples[,1:100,k], locs = which(val_data$location %in% LAD_NE), Replicates = rep_data$data)
k <- 11
rep_data <- ProcessData(subset(output_LAD, week == 12 & output == val_data$run$output[k]))
plot4 <- PlotReplicates(val_response_Het_rep$samples[,1:100,k], locs = which(val_data$location %in% LAD_NE), Replicates = rep_data$data)
cowplot::plot_grid(plot1,plot2,plot3,plot4)



#### Ward level ####
# Instead, attempt to model ward-level data, then compare aggregations
head(output_ward)
dim(output_ward) # 8071 wards x 700 simulations

# Use the same train/val split as above for comparison
n <- nrow(design)
train_prop <- 0.8
val_prop <- 1 - train_prop
set.seed(34810)
train_inds <- sample(1:n, train_prop*n)
train_labels <- design$output[train_inds]
train_data <- subset(output_ward, replicate == 1 & output %in% train_labels)
val_data <- subset(output_ward, replicate == 1 & !(output %in% train_labels))

train_data <- ProcessData(train_data, by_id = 'ward')
val_data <- ProcessData(val_data, by_id = 'ward')
dim(train_data$data) # 8071 x 200
dim(val_data$data) # 8071 x 50

# Construct basis
# basis_wards <- CountBasis(train_data$data, rank = 10) # rank is a design choice, use small value here for speed
# saveRDS(basis_wards, file = 'data/basis_wards_example.rds')
basis_wards <- readRDS('data/basis_wards_example.rds')
dim(basis_wards$tBasis) # 8071x10

# Plot leading basis vectors
PlotBasis(basis_wards, q = 1, LonLatWard[match(train_data$location, LonLatWard$ward),2:3])
PlotBasis(basis_wards, q = 2, LonLatWard[match(train_data$location, LonLatWard$ward),2:3])

# Combine inputs with coefficients
tData_ward <- GetEmDataCount(design_em[match(train_data$run$output, design_em$output),1:15], basis_wards, q = 10)
PlotPair(tData_ward, 'C1', 'C2', col = 'R0')
PlotPair(tData_ward, 'R0', 'ns', col = 'C1')

# Validation inputs
val_inputs <- subset(design_em, !(output %in% train_data$run$output))

# Emulate with hetGP
EmWardsHet <- BasisEmulatorsHet(tData_ward, 3, training_prop = 1)
Preds_val_Het <- BasisPredHet(val_inputs[,1:15], EmWardsHet)
val_wards_Het <- CountBasisEmSamples(Preds_val_Het, basis_wards, ReturnAll = TRUE)
ValidateSum(val_wards_Het$samples, val_data$data)

saveRDS(val_wards_Het, 'data/val_het_ward.rds')
saveRDS(val_data, 'data/val_data_ward.rds')

# Compare to LAD emulator
cowplot::plot_grid(ValidateSum(val_response_Het$samples, val_data$data), 
                   ValidateSum(val_wards_Het$samples, val_data$data), nrow = 1)
Metrics::rmse(apply(val_data$data,2,sum), apply(val_response$samples,3,sum))
Metrics::rmse(apply(val_data$data,2,sum), apply(val_response_Het$samples,3,sum))
Metrics::rmse(apply(val_data$data,2,sum), apply(val_wards_Het$samples,3,sum))

Metrics::rmse(log(apply(val_data$data,2,sum)+1), log(apply(val_response$samples,3,sum))+1)
Metrics::rmse(log(apply(val_data$data,2,sum)+1), log(apply(val_response_Het$samples,3,sum))+1)
Metrics::rmse(log(apply(val_data$data,2,sum)+1), log(apply(val_wards_Het$samples,3,sum))+1)

# NE samples - now wards rather than LADs
PlotSamplesCount(val_wards_Het$samples[,1:100,], locs = which(val_data$location %in% ward_NE), runs = 1:9, Truth = val_data$data)

# Generalises to replicates in the same way as above

