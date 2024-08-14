setwd("~/Dropbox/Exeter/Projects/UQ4Covid/CountBasis_paper")
source('code/1_ProcessData.R') # also loads in 0_CountBasis.R

# Plots for particular training/validation sets usually consistent with 2_WorkedExamples

#### Correlation between regional, overall totals ####
# North East is least correlated (on original and log scales), hence when plotting a particular region, will focus on this
output_total <- aggregate(cbind(cumH, deaths) ~ output + replicate, data = subset(output_LAD, week == 12), sum)
region_total <- aggregate(cbind(cumH, deaths) ~ output + replicate + region, data = subset(output_LAD, week == 12), sum)
regs <- unique(region_total$region)
for (rr in regs){
  print(c(rr, cor(output_total$deaths, subset(region_total, region == rr)$deaths)))
}
for (rr in regs){
  print(c(rr, cor(log(output_total$deaths+1), log(subset(region_total, region == rr)$deaths+1))))
}

#### Correlation between true data, emulator samples ####
# Not the best measure, particularly for lower values - output is non-linear
# val_data_LAD <- readRDS("data/val_data_LAD.rds") # validation data
# val_het_LAD <- readRDS("data/val_het_LAD.rds") # emulator samples for validation set
# lads <- val_data_LAD$location

# Select 2 LADs
# i <- 50; j <- 100
# total1 <- subset(output_LAD, week == 12 & LAD19CD == lads[i])
# total2 <- subset(output_LAD, week == 12 & LAD19CD == lads[j])
# cor(total1$deaths, total2$deaths)
# cor(c(val_het_LAD$samples[i,,]), c(val_het_LAD$samples[j,,]))
# cor(log(total1$deaths+1), log(total2$deaths+1))
# cor(c(log(val_het_LAD$samples[i,,]+1)), c(log(val_het_LAD$samples[j,,]+1)))
# 
# plot(total1$deaths, total2$deaths)
# plot(c(val_het_LAD$samples[i,,]), c(val_het_LAD$samples[j,,]))
# plot(c(log(val_het_LAD$samples[i,,]+1)), c(log(val_het_LAD$samples[j,,]+1)))

# Given the latent basis, we can write down correlations between f_i, f_j without needing to sample
# i.e. corr[f_i, f_j] = cov[f_i, f_j] / sd[f_i]*sd[f_j]
# For a particular input x
basis_deaths <- readRDS('data/basis_deaths_example.rds')
val_data_LAD <- readRDS("data/val_data_LAD.rds") # validation data
val_het_LAD <- readRDS("data/val_het_LAD.rds") # emulator samples for validation set
preds_het_LAD <- readRDS("data/preds_het_LAD.rds") # predictions on latent basis for validation set

# k <- 1
# z_mean <- as.numeric(basis_deaths$EnsembleMean + basis_deaths$LatentMean + basis_deaths$tBasis[,1:3] %*% preds_het_LAD$Expectation[k,])
# z_cov <- basis_deaths$tBasis[,1:3] %*% diag(preds_het_LAD$Variance[k,]) %*% t(basis_deaths$tBasis[,1:3])
# samps <- rmvnorm(10000, z_mean, z_cov) # 10000x339
# samps_p <- matrix(rpois(10000*339, exp(samps)), 10000, 339)
# 
# cor(samps_p[,1], samps_p[,300]) # sample correlation
# mm <- as.numeric(exp(z_mean + diag(z_cov)/2))
# vv <- mm + mm^2 * (exp(diag(z_cov)) - 1)
# cc <- mm[1]*mm*(exp(z_cov[1,])-1)
# cc[300] / sqrt(vv[1]*vv[300]) # theoretical correlation

# Select a LAD from each of the 9 regions
regions <- c('North East', 'North West', 'Yorkshire',
             'East Midlands', 'West Midlands', 'East of England',
             'London', 'South East', 'South West')
set.seed(5491)
k <- numeric(9)
k[1] <- sample(which(val_data_LAD$location %in% LAD_NE),1)
k[2] <- sample(which(val_data_LAD$location %in% LAD_NW),1)
k[3] <- sample(which(val_data_LAD$location %in% LAD_YO),1)
k[4] <- sample(which(val_data_LAD$location %in% LAD_EM),1)
k[5] <- sample(which(val_data_LAD$location %in% LAD_WM),1)
k[6] <- sample(which(val_data_LAD$location %in% LAD_EE),1)
k[7] <- sample(which(val_data_LAD$location %in% LAD_LO),1)
k[8] <- sample(which(val_data_LAD$location %in% LAD_SE),1)
k[9] <- sample(which(val_data_LAD$location %in% LAD_SW),1)
plot_data <- NULL
for (j in 1:9){
  ctrain <- cval <- csample <- csample2 <- numeric(339)
  for (i in 1:339){
    ctrain[i] <- cor(basis_deaths$Data[k[j],], basis_deaths$Data[i,])
    cval[i] <- cor(val_data_LAD$data[k[j],], val_data_LAD$data[i,])
    csample[i] <- cor(val_het_LAD$mean[k[j],], val_het_LAD$mean[i,])
    csample2[i] <- cor(c(val_het_LAD$samples[k[j],,]), c(val_het_LAD$samples[i,,]))
  }
  inds <- order(-ctrain)
  plot_data <- rbind(plot_data, data.frame(k = regions[j],
                                           LAD = 1:339,
                                           Training = ctrain[inds],
                                           Emulator = csample2[inds]))
                                           #Val = cval[inds],
                                           #Sample = csample[inds],
}
plot_data$k <- factor(plot_data$k, levels = regions)
plot_data <- melt(plot_data, id.vars = c('LAD', 'k'))

png('plots/correlations_by_region.png', width = 1344, height = 960, res = 200)
ggplot(plot_data, aes(LAD, value, col = variable)) + 
  geom_line(size = 0.5) +
  geom_line(data = subset(plot_data, variable == 'Training'), size = 0.75) + 
  facet_wrap(vars(k)) +
  scale_colour_manual(values = viridis(100)[c(71,21)]) +
  theme_bw() +
  labs(y = 'Correlation', col = 'Type')
dev.off()




#### Aggregated validation plots ####
val_data_LAD <- readRDS("data/val_data_LAD.rds") # validation data
val_het_LAD <- readRDS("data/val_het_LAD.rds") # emulator samples for validation set

plot1 <- ValidateSum(val_het_LAD$samples, val_data_LAD$data) + labs(title = '')
plot2 <- ValidateSum(val_het_LAD$samples, val_data_LAD$data, 
            locs = which(val_data_LAD$location %in% LAD_NE)) + labs(title = '')

val_data_ward <- readRDS("data/val_data_ward.rds") # validation data
val_het_ward <- readRDS("data/val_het_ward.rds") # emulator samples for validation set

plot3 <- ValidateSum(val_het_ward$samples, val_data_ward$data) + labs(title = '')
plot4 <- ValidateSum(val_het_ward$samples, val_data_ward$data, 
                     locs = which(val_data_ward$location %in% ward_NE)) + labs(title = '')

png('plots/validation_LAD_ward.png', width = 1344, height = 960, res = 200)
cowplot::plot_grid(plot1, plot3, plot2, plot4)
dev.off()

# These look very similar
# Might be more informative to directly plot LAD vs ward
LAD_total <- apply(val_het_LAD$mean, 2, sum)
ward_total <- apply(val_het_ward$mean, 2, sum)
plot_data <- data.frame(LAD = LAD_total, Ward = ward_total)
plot5 <- ggplot(plot_data, aes(Ward, LAD)) +
  geom_point() +
  geom_abline(slope = 1, alpha = 0.6) +
  scale_x_log10(breaks = c(0.1,1,10,100,1000,10000,100000,1000000), labels = c('0','1','10','100','1000','10000','100000','')) +
  scale_y_log10(breaks = c(0.1,1,10,100,1000,10000,100000,1000000), labels = c('0','1','10','100','1000','10000','100000','1000000')) +
  theme_bw() +
  theme(legend.position = 'none') +
  labs(x = 'PLNPCA-ward', y = 'PLNPCA-LAD')

# NE
LAD_total <- AggregateSamples(val_het_LAD$samples, locs = which(val_data_LAD$location %in% LAD_NE))
LAD_total <- apply(LAD_total$samples, 2, mean)
ward_total <- AggregateSamples(val_het_ward$samples, locs = which(val_data_ward$location %in% ward_NE))
ward_total <- apply(ward_total$samples, 2, mean)
plot_data <- data.frame(LAD = LAD_total, Ward = ward_total)
plot6 <- ggplot(plot_data, aes(Ward, LAD)) +
  geom_point() +
  geom_abline(slope = 1, alpha = 0.6) +
  scale_x_log10(breaks = c(0.1,1,10,100,1000,10000,100000,1000000), labels = c('0','1','10','100','1000','10000','100000','1000000')) +
  scale_y_log10(breaks = c(0.1,1,10,100,1000,10000,100000,1000000), labels = c('0','1','10','100','1000','10000','100000','1000000')) +
  theme_bw() +
  theme(legend.position = 'none') +
  labs(x = 'PLNPCA-ward', y = 'PLNPCA-LAD')

png('plots/validation_LAD_ward_v2.png', width = 1344, height = 960, res = 200)
cowplot::plot_grid(plot1, plot2, plot5, plot6, rel_heights = c(1,0.85))
dev.off()


# Repeating for other regions
val_data_LAD <- readRDS("data/val_data_LAD.rds") # validation data
val_het_LAD <- readRDS("data/val_het_LAD.rds") # emulator samples for validation set
val_data_ward <- readRDS("data/val_data_ward.rds") # validation data
val_het_ward <- readRDS("data/val_het_ward.rds") # emulator samples for validation set
plot1 <- ValidateSum(val_het_LAD$samples, val_data_LAD$data, 
                     locs = which(val_data_LAD$location %in% LAD_NW)) + labs(title = '')
plot2 <- ValidateSum(val_het_LAD$samples, val_data_LAD$data, 
                     locs = which(val_data_LAD$location %in% LAD_YO)) + labs(title = '')
LAD_total <- AggregateSamples(val_het_LAD$samples, locs = which(val_data_LAD$location %in% LAD_NW))
LAD_total <- apply(LAD_total$samples, 2, mean)
ward_total <- AggregateSamples(val_het_ward$samples, locs = which(val_data_ward$location %in% ward_NW))
ward_total <- apply(ward_total$samples, 2, mean)
plot_data <- data.frame(LAD = LAD_total, Ward = ward_total)
plot3 <- ggplot(plot_data, aes(Ward, LAD)) +
  geom_point() +
  geom_abline(slope = 1, alpha = 0.6) +
  scale_x_log10(breaks = c(0.1,1,10,100,1000,10000,100000,1000000), labels = c('0','1','10','100','1000','10000','','')) +
  scale_y_log10(breaks = c(0.1,1,10,100,1000,10000,100000,1000000), labels = c('0','1','10','100','1000','10000','100000','1000000')) +
  theme_bw() +
  theme(legend.position = 'none') +
  labs(x = 'PLNPCA-ward', y = 'PLNPCA-LAD')
LAD_total <- AggregateSamples(val_het_LAD$samples, locs = which(val_data_LAD$location %in% LAD_YO))
LAD_total <- apply(LAD_total$samples, 2, mean)
ward_total <- AggregateSamples(val_het_ward$samples, locs = which(val_data_ward$location %in% ward_YO))
ward_total <- apply(ward_total$samples, 2, mean)
plot_data <- data.frame(LAD = LAD_total, Ward = ward_total)
plot4 <- ggplot(plot_data, aes(Ward, LAD)) +
  geom_point() +
  geom_abline(slope = 1, alpha = 0.6) +
  scale_x_log10(breaks = c(0.1,1,10,100,1000,10000,100000,1000000), labels = c('0','1','10','100','1000','10000','100000','1000000')) +
  scale_y_log10(breaks = c(0.1,1,10,100,1000,10000,100000,1000000), labels = c('0','1','10','100','1000','10000','100000','1000000')) +
  theme_bw() +
  theme(legend.position = 'none') +
  labs(x = 'PLNPCA-ward', y = 'PLNPCA-LAD')

png('plots/validation_LAD_ward_NW_YO.png', width = 1344, height = 960, res = 200)
cowplot::plot_grid(plot1, plot2, plot3, plot4, rel_heights = c(1,0.85))
dev.off()

plot1 <- ValidateSum(val_het_LAD$samples, val_data_LAD$data, 
                     locs = which(val_data_LAD$location %in% LAD_EM)) + labs(title = '')
plot2 <- ValidateSum(val_het_LAD$samples, val_data_LAD$data, 
                     locs = which(val_data_LAD$location %in% LAD_WM)) + labs(title = '')
LAD_total <- AggregateSamples(val_het_LAD$samples, locs = which(val_data_LAD$location %in% LAD_EM))
LAD_total <- apply(LAD_total$samples, 2, mean)
ward_total <- AggregateSamples(val_het_ward$samples, locs = which(val_data_ward$location %in% ward_EM))
ward_total <- apply(ward_total$samples, 2, mean)
plot_data <- data.frame(LAD = LAD_total, Ward = ward_total)
plot3 <- ggplot(plot_data, aes(Ward, LAD)) +
  geom_point() +
  geom_abline(slope = 1, alpha = 0.6) +
  scale_x_log10(breaks = c(0.1,1,10,100,1000,10000,100000,1000000), labels = c('0','1','10','100','1000','10000','','')) +
  scale_y_log10(breaks = c(0.1,1,10,100,1000,10000,100000,1000000), labels = c('0','1','10','100','1000','10000','100000','1000000')) +
  theme_bw() +
  theme(legend.position = 'none') +
  labs(x = 'PLNPCA-ward', y = 'PLNPCA-LAD')
LAD_total <- AggregateSamples(val_het_LAD$samples, locs = which(val_data_LAD$location %in% LAD_WM))
LAD_total <- apply(LAD_total$samples, 2, mean)
ward_total <- AggregateSamples(val_het_ward$samples, locs = which(val_data_ward$location %in% ward_WM))
ward_total <- apply(ward_total$samples, 2, mean)
plot_data <- data.frame(LAD = LAD_total, Ward = ward_total)
plot4 <- ggplot(plot_data, aes(Ward, LAD)) +
  geom_point() +
  geom_abline(slope = 1, alpha = 0.6) +
  scale_x_log10(breaks = c(0.1,1,10,100,1000,10000,100000,1000000), labels = c('0','1','10','100','1000','10000','','1000000')) +
  scale_y_log10(breaks = c(0.1,1,10,100,1000,10000,100000,1000000), labels = c('0','1','10','100','1000','10000','100000','1000000')) +
  theme_bw() +
  theme(legend.position = 'none') +
  labs(x = 'PLNPCA-ward', y = 'PLNPCA-LAD')

png('plots/validation_LAD_ward_EM_WM.png', width = 1344, height = 960, res = 200)
cowplot::plot_grid(plot1, plot2, plot3, plot4, rel_heights = c(1,0.85))
dev.off()


plot1 <- ValidateSum(val_het_LAD$samples, val_data_LAD$data, 
                     locs = which(val_data_LAD$location %in% LAD_EE)) + labs(title = '')
plot2 <- ValidateSum(val_het_LAD$samples, val_data_LAD$data, 
                     locs = which(val_data_LAD$location %in% LAD_LO)) + labs(title = '')
LAD_total <- AggregateSamples(val_het_LAD$samples, locs = which(val_data_LAD$location %in% LAD_EE))
LAD_total <- apply(LAD_total$samples, 2, mean)
ward_total <- AggregateSamples(val_het_ward$samples, locs = which(val_data_ward$location %in% ward_EE))
ward_total <- apply(ward_total$samples, 2, mean)
plot_data <- data.frame(LAD = LAD_total, Ward = ward_total)
plot3 <- ggplot(plot_data, aes(Ward, LAD)) +
  geom_point() +
  geom_abline(slope = 1, alpha = 0.6) +
  scale_x_log10(breaks = c(0.1,1,10,100,1000,10000,100000,1000000), labels = c('0','1','10','100','1000','10000','100000','')) +
  scale_y_log10(breaks = c(0.1,1,10,100,1000,10000,100000,1000000), labels = c('0','1','10','100','1000','10000','100000','1000000')) +
  theme_bw() +
  theme(legend.position = 'none') +
  labs(x = 'PLNPCA-ward', y = 'PLNPCA-LAD')
LAD_total <- AggregateSamples(val_het_LAD$samples, locs = which(val_data_LAD$location %in% LAD_LO))
LAD_total <- apply(LAD_total$samples, 2, mean)
ward_total <- AggregateSamples(val_het_ward$samples, locs = which(val_data_ward$location %in% ward_LO))
ward_total <- apply(ward_total$samples, 2, mean)
plot_data <- data.frame(LAD = LAD_total, Ward = ward_total)
plot4 <- ggplot(plot_data, aes(Ward, LAD)) +
  geom_point() +
  geom_abline(slope = 1, alpha = 0.6) +
  scale_x_log10(breaks = c(0.1,1,10,100,1000,10000,100000,1000000), labels = c('0','1','10','100','1000','10000','100000','1000000')) +
  scale_y_log10(breaks = c(0.1,1,10,100,1000,10000,100000,1000000), labels = c('0','1','10','100','1000','10000','100000','1000000')) +
  theme_bw() +
  theme(legend.position = 'none') +
  labs(x = 'PLNPCA-ward', y = 'PLNPCA-LAD')

png('plots/validation_LAD_ward_EE_LO.png', width = 1344, height = 960, res = 200)
cowplot::plot_grid(plot1, plot2, plot3, plot4, rel_heights = c(1,0.85))
dev.off()


plot1 <- ValidateSum(val_het_LAD$samples, val_data_LAD$data, 
                     locs = which(val_data_LAD$location %in% LAD_SE)) + labs(title = '')
plot2 <- ValidateSum(val_het_LAD$samples, val_data_LAD$data, 
                     locs = which(val_data_LAD$location %in% LAD_SW)) + labs(title = '')
LAD_total <- AggregateSamples(val_het_LAD$samples, locs = which(val_data_LAD$location %in% LAD_SE))
LAD_total <- apply(LAD_total$samples, 2, mean)
ward_total <- AggregateSamples(val_het_ward$samples, locs = which(val_data_ward$location %in% ward_SE))
ward_total <- apply(ward_total$samples, 2, mean)
plot_data <- data.frame(LAD = LAD_total, Ward = ward_total)
plot3 <- ggplot(plot_data, aes(Ward, LAD)) +
  geom_point() +
  geom_abline(slope = 1, alpha = 0.6) +
  scale_x_log10(breaks = c(0.1,1,10,100,1000,10000,100000,1000000), labels = c('0','1','10','100','1000','10000','100000','')) +
  scale_y_log10(breaks = c(0.1,1,10,100,1000,10000,100000,1000000), labels = c('0','1','10','100','1000','10000','100000','1000000')) +
  theme_bw() +
  theme(legend.position = 'none') +
  labs(x = 'PLNPCA-ward', y = 'PLNPCA-LAD')
LAD_total <- AggregateSamples(val_het_LAD$samples, locs = which(val_data_LAD$location %in% LAD_SW))
LAD_total <- apply(LAD_total$samples, 2, mean)
ward_total <- AggregateSamples(val_het_ward$samples, locs = which(val_data_ward$location %in% ward_SW))
ward_total <- apply(ward_total$samples, 2, mean)
plot_data <- data.frame(LAD = LAD_total, Ward = ward_total)
plot4 <- ggplot(plot_data, aes(Ward, LAD)) +
  geom_point() +
  geom_abline(slope = 1, alpha = 0.6) +
  scale_x_log10(breaks = c(0.1,1,10,100,1000,10000,100000,1000000), labels = c('0','1','10','100','1000','10000','','1000000')) +
  scale_y_log10(breaks = c(0.1,1,10,100,1000,10000,100000,1000000), labels = c('0','1','10','100','1000','10000','100000','1000000')) +
  theme_bw() +
  theme(legend.position = 'none') +
  labs(x = 'PLNPCA-ward', y = 'PLNPCA-LAD')

png('plots/validation_LAD_ward_SE_SW.png', width = 1344, height = 960, res = 200)
cowplot::plot_grid(plot1, plot2, plot3, plot4, rel_heights = c(1,0.85))
dev.off()




#### With reps ####
CoeffPlots <- readRDS('data/CoeffPlots.rds')
png('plots/coeffs_reps.png', width = 1344, height = 960, res = 200)
cowplot::plot_grid(CoeffPlots[[1]],CoeffPlots[[2]],CoeffPlots[[3]],CoeffPlots[[4]])
dev.off()

val_data_LAD_reps <- readRDS("data/val_data_LAD_reps.rds") # validation data
val_het_LAD_reps <- readRDS("data/val_het_LAD_reps.rds") # emulator samples for validation set

plot1 <- ValidateSum(val_het_LAD_reps$samples, val_data_LAD_reps$data) + labs(title = '')
plot2 <- ValidateSum(val_het_LAD_reps$samples, val_data_LAD_reps$data, 
                     locs = which(val_data_LAD_reps$location %in% LAD_NE)) + labs(title = '')

png('plots/validation_LAD_reps.png', width = 1344, height = 600, res = 200)
cowplot::plot_grid(plot1, plot2)
dev.off()

plot1 <- ValidateSum(val_het_LAD_reps$samples, val_data_LAD_reps$data, 
                     locs = which(val_data_LAD_reps$location %in% LAD_NW)) + labs(title = '')
plot2 <- ValidateSum(val_het_LAD_reps$samples, val_data_LAD_reps$data, 
                     locs = which(val_data_LAD_reps$location %in% LAD_YO)) + labs(title = '')
png('plots/validation_LAD_reps_NW_YO.png', width = 1344, height = 600, res = 200)
cowplot::plot_grid(plot1, plot2)
dev.off()

plot1 <- ValidateSum(val_het_LAD_reps$samples, val_data_LAD_reps$data, 
                     locs = which(val_data_LAD_reps$location %in% LAD_EM)) + labs(title = '')
plot2 <- ValidateSum(val_het_LAD_reps$samples, val_data_LAD_reps$data, 
                     locs = which(val_data_LAD_reps$location %in% LAD_WM)) + labs(title = '')
png('plots/validation_LAD_reps_EM_WM.png', width = 1344, height = 600, res = 200)
cowplot::plot_grid(plot1, plot2)
dev.off()

plot1 <- ValidateSum(val_het_LAD_reps$samples, val_data_LAD_reps$data, 
                     locs = which(val_data_LAD_reps$location %in% LAD_EE)) + labs(title = '')
plot2 <- ValidateSum(val_het_LAD_reps$samples, val_data_LAD_reps$data, 
                     locs = which(val_data_LAD_reps$location %in% LAD_LO)) + labs(title = '') +
  scale_x_log10(breaks = c(0.1,1,10,100,1000,10000,100000,1000000), labels = c('0','1','10','100','1000','10000','','1000000')) 

png('plots/validation_LAD_reps_EE_LO.png', width = 1344, height = 600, res = 200)
cowplot::plot_grid(plot1, plot2)
dev.off()

plot1 <- ValidateSum(val_het_LAD_reps$samples, val_data_LAD_reps$data, 
                     locs = which(val_data_LAD_reps$location %in% LAD_SE)) + labs(title = '') +
  scale_x_log10(breaks = c(0.1,1,10,100,1000,10000,100000,1000000), labels = c('0','1','10','100','1000','10000','','1000000')) 
plot2 <- ValidateSum(val_het_LAD_reps$samples, val_data_LAD_reps$data, 
                     locs = which(val_data_LAD_reps$location %in% LAD_SW)) + labs(title = '')
png('plots/validation_LAD_reps_SE_SW.png', width = 1344, height = 600, res = 200)
cowplot::plot_grid(plot1, plot2)
dev.off()



#### Single location validation ####
ValidateSum(val_het_LAD$samples, val_data_LAD$data, locs = 49) + labs(title = '')




#### Plotting samples ####
png('plots/samples_het_NE.png', width = 1344, height = 960, res = 200)
PlotSamplesCount(val_het_LAD$samples, locs = which(val_data_LAD$location %in% LAD_NE), runs = 1:9, Truth = val_data_LAD$data) + 
  labs(x = 'LAD', y = 'Deaths')
dev.off()

png('plots/samples_het_NW.png', width = 1344, height = 960, res = 200)
PlotSamplesCount(val_het_LAD$samples, locs = which(val_data_LAD$location %in% LAD_NW), runs = 1:9, Truth = val_data_LAD$data) + 
  labs(x = 'LAD', y = 'Deaths')
dev.off()

png('plots/samples_het_YO.png', width = 1344, height = 960, res = 200)
PlotSamplesCount(val_het_LAD$samples, locs = which(val_data_LAD$location %in% LAD_YO), runs = 1:9, Truth = val_data_LAD$data) + 
  labs(x = 'LAD', y = 'Deaths')
dev.off()

png('plots/samples_het_EM.png', width = 1344, height = 960, res = 200)
PlotSamplesCount(val_het_LAD$samples, locs = which(val_data_LAD$location %in% LAD_EM), runs = 1:9, Truth = val_data_LAD$data) + 
  labs(x = 'LAD', y = 'Deaths')
dev.off()

png('plots/samples_het_WM.png', width = 1344, height = 960, res = 200)
PlotSamplesCount(val_het_LAD$samples, locs = which(val_data_LAD$location %in% LAD_WM), runs = 1:9, Truth = val_data_LAD$data) + 
  labs(x = 'LAD', y = 'Deaths')
dev.off()

png('plots/samples_het_EE.png', width = 1344, height = 960, res = 200)
PlotSamplesCount(val_het_LAD$samples, locs = which(val_data_LAD$location %in% LAD_EE), runs = 1:9, Truth = val_data_LAD$data) + 
  labs(x = 'LAD', y = 'Deaths')
dev.off()

png('plots/samples_het_LO.png', width = 1344, height = 960, res = 200)
PlotSamplesCount(val_het_LAD$samples, locs = which(val_data_LAD$location %in% LAD_LO), runs = 1:9, Truth = val_data_LAD$data) + 
  labs(x = 'LAD', y = 'Deaths')
dev.off()

png('plots/samples_het_SE.png', width = 1344, height = 960, res = 200)
PlotSamplesCount(val_het_LAD$samples, locs = which(val_data_LAD$location %in% LAD_SE), runs = 1:9, Truth = val_data_LAD$data) + 
  labs(x = 'LAD', y = 'Deaths')
dev.off()

png('plots/samples_het_SW.png', width = 1344, height = 960, res = 200)
PlotSamplesCount(val_het_LAD$samples, locs = which(val_data_LAD$location %in% LAD_SW), runs = 1:9, Truth = val_data_LAD$data) + 
  labs(x = 'LAD', y = 'Deaths')
dev.off()


# Ward-level
png('plots/samples_het_NE_ward.png', width = 1344, height = 960, res = 200)
PlotSamplesCount(val_het_ward$samples[,1:100,], locs = which(val_data_ward$location %in% ward_NE)[1:100], runs = 1:9, Truth = val_data_ward$data) + 
  labs(x = 'Ward', y = 'Deaths')
dev.off()

png('plots/samples_het_NW_ward.png', width = 1344, height = 960, res = 200)
PlotSamplesCount(val_het_ward$samples[,1:100,], locs = which(val_data_ward$location %in% ward_NW)[1:100], runs = 1:9, Truth = val_data_ward$data) + 
  labs(x = 'Ward', y = 'Deaths')
dev.off()

png('plots/samples_het_YO_ward.png', width = 1344, height = 960, res = 200)
PlotSamplesCount(val_het_ward$samples[,1:100,], locs = which(val_data_ward$location %in% ward_YO)[1:100], runs = 1:9, Truth = val_data_ward$data) + 
  labs(x = 'Ward', y = 'Deaths')
dev.off()

png('plots/samples_het_EM_ward.png', width = 1344, height = 960, res = 200)
PlotSamplesCount(val_het_ward$samples[,1:100,], locs = which(val_data_ward$location %in% ward_EM)[1:100], runs = 1:9, Truth = val_data_ward$data) + 
  labs(x = 'Ward', y = 'Deaths')
dev.off()

png('plots/samples_het_WM_ward.png', width = 1344, height = 960, res = 200)
PlotSamplesCount(val_het_ward$samples[,1:100,], locs = which(val_data_ward$location %in% ward_WM)[1:100], runs = 1:9, Truth = val_data_ward$data) + 
  labs(x = 'Ward', y = 'Deaths')
dev.off()

png('plots/samples_het_EE_ward.png', width = 1344, height = 960, res = 200)
PlotSamplesCount(val_het_ward$samples[,1:100,], locs = which(val_data_ward$location %in% ward_EE)[1:100], runs = 1:9, Truth = val_data_ward$data) + 
  labs(x = 'Ward', y = 'Deaths')
dev.off()

png('plots/samples_het_LO_ward.png', width = 1344, height = 960, res = 200)
PlotSamplesCount(val_het_ward$samples[,1:100,], locs = which(val_data_ward$location %in% ward_LO)[1:100], runs = 1:9, Truth = val_data_ward$data) + 
  labs(x = 'Ward', y = 'Deaths')
dev.off()

png('plots/samples_het_SE_ward.png', width = 1344, height = 960, res = 200)
PlotSamplesCount(val_het_ward$samples[,1:100,], locs = which(val_data_ward$location %in% ward_SE)[1:100], runs = 1:9, Truth = val_data_ward$data) + 
  labs(x = 'Ward', y = 'Deaths')
dev.off()

png('plots/samples_het_SW_ward.png', width = 1344, height = 960, res = 200)
PlotSamplesCount(val_het_ward$samples[,1:100,], locs = which(val_data_ward$location %in% ward_SW)[1:100], runs = 1:9, Truth = val_data_ward$data) + 
  labs(x = 'Ward', y = 'Deaths')
dev.off()




# With replicates
png('plots/samples_het_NE_reps.png', width = 1344, height = 960, res = 200)
k <- 1;rep_data <- ProcessData(subset(output_LAD, week == 12 & output == val_data_LAD_reps$run$output[k]))
plot1 <- PlotReplicates(val_het_LAD_reps$samples[,1:100,k], locs = which(val_data_LAD_reps$location %in% LAD_NE), Replicates = rep_data$data) + 
  labs(x = 'LAD', y = 'Deaths')
k <- 2;rep_data <- ProcessData(subset(output_LAD, week == 12 & output == val_data_LAD_reps$run$output[k]))
plot2 <- PlotReplicates(val_het_LAD_reps$samples[,1:100,k], locs = which(val_data_LAD_reps$location %in% LAD_NE), Replicates = rep_data$data) + 
  labs(x = 'LAD', y = 'Deaths')
k <- 3;rep_data <- ProcessData(subset(output_LAD, week == 12 & output == val_data_LAD_reps$run$output[k]))
plot3 <- PlotReplicates(val_het_LAD_reps$samples[,1:100,k], locs = which(val_data_LAD_reps$location %in% LAD_NE), Replicates = rep_data$data) + 
  labs(x = 'LAD', y = 'Deaths')
k <- 4;rep_data <- ProcessData(subset(output_LAD, week == 12 & output == val_data_LAD_reps$run$output[k]))
plot4 <- PlotReplicates(val_het_LAD_reps$samples[,1:100,k], locs = which(val_data_LAD_reps$location %in% LAD_NE), Replicates = rep_data$data) + 
  labs(x = 'LAD', y = 'Deaths')
k <- 5;rep_data <- ProcessData(subset(output_LAD, week == 12 & output == val_data_LAD_reps$run$output[k]))
plot5 <- PlotReplicates(val_het_LAD_reps$samples[,1:100,k], locs = which(val_data_LAD_reps$location %in% LAD_NE), Replicates = rep_data$data) + 
  labs(x = 'LAD', y = 'Deaths')
k <- 6;rep_data <- ProcessData(subset(output_LAD, week == 12 & output == val_data_LAD_reps$run$output[k]))
plot6 <- PlotReplicates(val_het_LAD_reps$samples[,1:100,k], locs = which(val_data_LAD_reps$location %in% LAD_NE), Replicates = rep_data$data) + 
  labs(x = 'LAD', y = 'Deaths')
cowplot::plot_grid(plot1,plot2,plot3,plot4,plot5,plot6)
dev.off()

png('plots/samples_het_NW_reps.png', width = 1344, height = 960, res = 200)
k <- 1;rep_data <- ProcessData(subset(output_LAD, week == 12 & output == val_data_LAD_reps$run$output[k]))
plot1 <- PlotReplicates(val_het_LAD_reps$samples[,1:100,k], locs = which(val_data_LAD_reps$location %in% LAD_NW), Replicates = rep_data$data) + 
  labs(x = 'LAD', y = 'Deaths')
k <- 2;rep_data <- ProcessData(subset(output_LAD, week == 12 & output == val_data_LAD_reps$run$output[k]))
plot2 <- PlotReplicates(val_het_LAD_reps$samples[,1:100,k], locs = which(val_data_LAD_reps$location %in% LAD_NW), Replicates = rep_data$data) + 
  labs(x = 'LAD', y = 'Deaths')
k <- 3;rep_data <- ProcessData(subset(output_LAD, week == 12 & output == val_data_LAD_reps$run$output[k]))
plot3 <- PlotReplicates(val_het_LAD_reps$samples[,1:100,k], locs = which(val_data_LAD_reps$location %in% LAD_NW), Replicates = rep_data$data) + 
  labs(x = 'LAD', y = 'Deaths')
k <- 4;rep_data <- ProcessData(subset(output_LAD, week == 12 & output == val_data_LAD_reps$run$output[k]))
plot4 <- PlotReplicates(val_het_LAD_reps$samples[,1:100,k], locs = which(val_data_LAD_reps$location %in% LAD_NW), Replicates = rep_data$data) + 
  labs(x = 'LAD', y = 'Deaths')
k <- 5;rep_data <- ProcessData(subset(output_LAD, week == 12 & output == val_data_LAD_reps$run$output[k]))
plot5 <- PlotReplicates(val_het_LAD_reps$samples[,1:100,k], locs = which(val_data_LAD_reps$location %in% LAD_NW), Replicates = rep_data$data) + 
  labs(x = 'LAD', y = 'Deaths')
k <- 6;rep_data <- ProcessData(subset(output_LAD, week == 12 & output == val_data_LAD_reps$run$output[k]))
plot6 <- PlotReplicates(val_het_LAD_reps$samples[,1:100,k], locs = which(val_data_LAD_reps$location %in% LAD_NW), Replicates = rep_data$data) + 
  labs(x = 'LAD', y = 'Deaths')
cowplot::plot_grid(plot1,plot2,plot3,plot4,plot5,plot6)
dev.off()

png('plots/samples_het_YO_reps.png', width = 1344, height = 960, res = 200)
k <- 1;rep_data <- ProcessData(subset(output_LAD, week == 12 & output == val_data_LAD_reps$run$output[k]))
plot1 <- PlotReplicates(val_het_LAD_reps$samples[,1:100,k], locs = which(val_data_LAD_reps$location %in% LAD_YO), Replicates = rep_data$data) + 
  labs(x = 'LAD', y = 'Deaths')
k <- 2;rep_data <- ProcessData(subset(output_LAD, week == 12 & output == val_data_LAD_reps$run$output[k]))
plot2 <- PlotReplicates(val_het_LAD_reps$samples[,1:100,k], locs = which(val_data_LAD_reps$location %in% LAD_YO), Replicates = rep_data$data) + 
  labs(x = 'LAD', y = 'Deaths')
k <- 3;rep_data <- ProcessData(subset(output_LAD, week == 12 & output == val_data_LAD_reps$run$output[k]))
plot3 <- PlotReplicates(val_het_LAD_reps$samples[,1:100,k], locs = which(val_data_LAD_reps$location %in% LAD_YO), Replicates = rep_data$data) + 
  labs(x = 'LAD', y = 'Deaths')
k <- 4;rep_data <- ProcessData(subset(output_LAD, week == 12 & output == val_data_LAD_reps$run$output[k]))
plot4 <- PlotReplicates(val_het_LAD_reps$samples[,1:100,k], locs = which(val_data_LAD_reps$location %in% LAD_YO), Replicates = rep_data$data) + 
  labs(x = 'LAD', y = 'Deaths')
k <- 5;rep_data <- ProcessData(subset(output_LAD, week == 12 & output == val_data_LAD_reps$run$output[k]))
plot5 <- PlotReplicates(val_het_LAD_reps$samples[,1:100,k], locs = which(val_data_LAD_reps$location %in% LAD_YO), Replicates = rep_data$data) + 
  labs(x = 'LAD', y = 'Deaths')
k <- 6;rep_data <- ProcessData(subset(output_LAD, week == 12 & output == val_data_LAD_reps$run$output[k]))
plot6 <- PlotReplicates(val_het_LAD_reps$samples[,1:100,k], locs = which(val_data_LAD_reps$location %in% LAD_YO), Replicates = rep_data$data) + 
  labs(x = 'LAD', y = 'Deaths')
cowplot::plot_grid(plot1,plot2,plot3,plot4,plot5,plot6)
dev.off()

png('plots/samples_het_EM_reps.png', width = 1344, height = 960, res = 200)
k <- 1;rep_data <- ProcessData(subset(output_LAD, week == 12 & output == val_data_LAD_reps$run$output[k]))
plot1 <- PlotReplicates(val_het_LAD_reps$samples[,1:100,k], locs = which(val_data_LAD_reps$location %in% LAD_EM), Replicates = rep_data$data) + 
  labs(x = 'LAD', y = 'Deaths')
k <- 2;rep_data <- ProcessData(subset(output_LAD, week == 12 & output == val_data_LAD_reps$run$output[k]))
plot2 <- PlotReplicates(val_het_LAD_reps$samples[,1:100,k], locs = which(val_data_LAD_reps$location %in% LAD_EM), Replicates = rep_data$data) + 
  labs(x = 'LAD', y = 'Deaths')
k <- 3;rep_data <- ProcessData(subset(output_LAD, week == 12 & output == val_data_LAD_reps$run$output[k]))
plot3 <- PlotReplicates(val_het_LAD_reps$samples[,1:100,k], locs = which(val_data_LAD_reps$location %in% LAD_EM), Replicates = rep_data$data) + 
  labs(x = 'LAD', y = 'Deaths')
k <- 4;rep_data <- ProcessData(subset(output_LAD, week == 12 & output == val_data_LAD_reps$run$output[k]))
plot4 <- PlotReplicates(val_het_LAD_reps$samples[,1:100,k], locs = which(val_data_LAD_reps$location %in% LAD_EM), Replicates = rep_data$data) + 
  labs(x = 'LAD', y = 'Deaths')
k <- 5;rep_data <- ProcessData(subset(output_LAD, week == 12 & output == val_data_LAD_reps$run$output[k]))
plot5 <- PlotReplicates(val_het_LAD_reps$samples[,1:100,k], locs = which(val_data_LAD_reps$location %in% LAD_EM), Replicates = rep_data$data) + 
  labs(x = 'LAD', y = 'Deaths')
k <- 6;rep_data <- ProcessData(subset(output_LAD, week == 12 & output == val_data_LAD_reps$run$output[k]))
plot6 <- PlotReplicates(val_het_LAD_reps$samples[,1:100,k], locs = which(val_data_LAD_reps$location %in% LAD_EM), Replicates = rep_data$data) + 
  labs(x = 'LAD', y = 'Deaths')
cowplot::plot_grid(plot1,plot2,plot3,plot4,plot5,plot6)
dev.off()

png('plots/samples_het_WM_reps.png', width = 1344, height = 960, res = 200)
k <- 1;rep_data <- ProcessData(subset(output_LAD, week == 12 & output == val_data_LAD_reps$run$output[k]))
plot1 <- PlotReplicates(val_het_LAD_reps$samples[,1:100,k], locs = which(val_data_LAD_reps$location %in% LAD_WM), Replicates = rep_data$data) + 
  labs(x = 'LAD', y = 'Deaths')
k <- 2;rep_data <- ProcessData(subset(output_LAD, week == 12 & output == val_data_LAD_reps$run$output[k]))
plot2 <- PlotReplicates(val_het_LAD_reps$samples[,1:100,k], locs = which(val_data_LAD_reps$location %in% LAD_WM), Replicates = rep_data$data) + 
  labs(x = 'LAD', y = 'Deaths')
k <- 3;rep_data <- ProcessData(subset(output_LAD, week == 12 & output == val_data_LAD_reps$run$output[k]))
plot3 <- PlotReplicates(val_het_LAD_reps$samples[,1:100,k], locs = which(val_data_LAD_reps$location %in% LAD_WM), Replicates = rep_data$data) + 
  labs(x = 'LAD', y = 'Deaths')
k <- 4;rep_data <- ProcessData(subset(output_LAD, week == 12 & output == val_data_LAD_reps$run$output[k]))
plot4 <- PlotReplicates(val_het_LAD_reps$samples[,1:100,k], locs = which(val_data_LAD_reps$location %in% LAD_WM), Replicates = rep_data$data) + 
  labs(x = 'LAD', y = 'Deaths')
k <- 5;rep_data <- ProcessData(subset(output_LAD, week == 12 & output == val_data_LAD_reps$run$output[k]))
plot5 <- PlotReplicates(val_het_LAD_reps$samples[,1:100,k], locs = which(val_data_LAD_reps$location %in% LAD_WM), Replicates = rep_data$data) + 
  labs(x = 'LAD', y = 'Deaths')
k <- 6;rep_data <- ProcessData(subset(output_LAD, week == 12 & output == val_data_LAD_reps$run$output[k]))
plot6 <- PlotReplicates(val_het_LAD_reps$samples[,1:100,k], locs = which(val_data_LAD_reps$location %in% LAD_WM), Replicates = rep_data$data) + 
  labs(x = 'LAD', y = 'Deaths')
cowplot::plot_grid(plot1,plot2,plot3,plot4,plot5,plot6)
dev.off()

png('plots/samples_het_EE_reps.png', width = 1344, height = 960, res = 200)
k <- 1;rep_data <- ProcessData(subset(output_LAD, week == 12 & output == val_data_LAD_reps$run$output[k]))
plot1 <- PlotReplicates(val_het_LAD_reps$samples[,1:100,k], locs = which(val_data_LAD_reps$location %in% LAD_EE), Replicates = rep_data$data) + 
  labs(x = 'LAD', y = 'Deaths')
k <- 2;rep_data <- ProcessData(subset(output_LAD, week == 12 & output == val_data_LAD_reps$run$output[k]))
plot2 <- PlotReplicates(val_het_LAD_reps$samples[,1:100,k], locs = which(val_data_LAD_reps$location %in% LAD_EE), Replicates = rep_data$data) + 
  labs(x = 'LAD', y = 'Deaths')
k <- 3;rep_data <- ProcessData(subset(output_LAD, week == 12 & output == val_data_LAD_reps$run$output[k]))
plot3 <- PlotReplicates(val_het_LAD_reps$samples[,1:100,k], locs = which(val_data_LAD_reps$location %in% LAD_EE), Replicates = rep_data$data) + 
  labs(x = 'LAD', y = 'Deaths')
k <- 4;rep_data <- ProcessData(subset(output_LAD, week == 12 & output == val_data_LAD_reps$run$output[k]))
plot4 <- PlotReplicates(val_het_LAD_reps$samples[,1:100,k], locs = which(val_data_LAD_reps$location %in% LAD_EE), Replicates = rep_data$data) + 
  labs(x = 'LAD', y = 'Deaths')
k <- 5;rep_data <- ProcessData(subset(output_LAD, week == 12 & output == val_data_LAD_reps$run$output[k]))
plot5 <- PlotReplicates(val_het_LAD_reps$samples[,1:100,k], locs = which(val_data_LAD_reps$location %in% LAD_EE), Replicates = rep_data$data) + 
  labs(x = 'LAD', y = 'Deaths')
k <- 6;rep_data <- ProcessData(subset(output_LAD, week == 12 & output == val_data_LAD_reps$run$output[k]))
plot6 <- PlotReplicates(val_het_LAD_reps$samples[,1:100,k], locs = which(val_data_LAD_reps$location %in% LAD_EE), Replicates = rep_data$data) + 
  labs(x = 'LAD', y = 'Deaths')
cowplot::plot_grid(plot1,plot2,plot3,plot4,plot5,plot6)
dev.off()

png('plots/samples_het_LO_reps.png', width = 1344, height = 960, res = 200)
k <- 1;rep_data <- ProcessData(subset(output_LAD, week == 12 & output == val_data_LAD_reps$run$output[k]))
plot1 <- PlotReplicates(val_het_LAD_reps$samples[,1:100,k], locs = which(val_data_LAD_reps$location %in% LAD_LO), Replicates = rep_data$data) + 
  labs(x = 'LAD', y = 'Deaths')
k <- 2;rep_data <- ProcessData(subset(output_LAD, week == 12 & output == val_data_LAD_reps$run$output[k]))
plot2 <- PlotReplicates(val_het_LAD_reps$samples[,1:100,k], locs = which(val_data_LAD_reps$location %in% LAD_LO), Replicates = rep_data$data) + 
  labs(x = 'LAD', y = 'Deaths')
k <- 3;rep_data <- ProcessData(subset(output_LAD, week == 12 & output == val_data_LAD_reps$run$output[k]))
plot3 <- PlotReplicates(val_het_LAD_reps$samples[,1:100,k], locs = which(val_data_LAD_reps$location %in% LAD_LO), Replicates = rep_data$data) + 
  labs(x = 'LAD', y = 'Deaths')
k <- 4;rep_data <- ProcessData(subset(output_LAD, week == 12 & output == val_data_LAD_reps$run$output[k]))
plot4 <- PlotReplicates(val_het_LAD_reps$samples[,1:100,k], locs = which(val_data_LAD_reps$location %in% LAD_LO), Replicates = rep_data$data) + 
  labs(x = 'LAD', y = 'Deaths')
k <- 5;rep_data <- ProcessData(subset(output_LAD, week == 12 & output == val_data_LAD_reps$run$output[k]))
plot5 <- PlotReplicates(val_het_LAD_reps$samples[,1:100,k], locs = which(val_data_LAD_reps$location %in% LAD_LO), Replicates = rep_data$data) + 
  labs(x = 'LAD', y = 'Deaths')
k <- 6;rep_data <- ProcessData(subset(output_LAD, week == 12 & output == val_data_LAD_reps$run$output[k]))
plot6 <- PlotReplicates(val_het_LAD_reps$samples[,1:100,k], locs = which(val_data_LAD_reps$location %in% LAD_LO), Replicates = rep_data$data) + 
  labs(x = 'LAD', y = 'Deaths')
cowplot::plot_grid(plot1,plot2,plot3,plot4,plot5,plot6)
dev.off()

png('plots/samples_het_SE_reps.png', width = 1344, height = 960, res = 200)
k <- 1;rep_data <- ProcessData(subset(output_LAD, week == 12 & output == val_data_LAD_reps$run$output[k]))
plot1 <- PlotReplicates(val_het_LAD_reps$samples[,1:100,k], locs = which(val_data_LAD_reps$location %in% LAD_SE), Replicates = rep_data$data) + 
  labs(x = 'LAD', y = 'Deaths')
k <- 2;rep_data <- ProcessData(subset(output_LAD, week == 12 & output == val_data_LAD_reps$run$output[k]))
plot2 <- PlotReplicates(val_het_LAD_reps$samples[,1:100,k], locs = which(val_data_LAD_reps$location %in% LAD_SE), Replicates = rep_data$data) + 
  labs(x = 'LAD', y = 'Deaths')
k <- 3;rep_data <- ProcessData(subset(output_LAD, week == 12 & output == val_data_LAD_reps$run$output[k]))
plot3 <- PlotReplicates(val_het_LAD_reps$samples[,1:100,k], locs = which(val_data_LAD_reps$location %in% LAD_SE), Replicates = rep_data$data) + 
  labs(x = 'LAD', y = 'Deaths')
k <- 4;rep_data <- ProcessData(subset(output_LAD, week == 12 & output == val_data_LAD_reps$run$output[k]))
plot4 <- PlotReplicates(val_het_LAD_reps$samples[,1:100,k], locs = which(val_data_LAD_reps$location %in% LAD_SE), Replicates = rep_data$data) + 
  labs(x = 'LAD', y = 'Deaths')
k <- 5;rep_data <- ProcessData(subset(output_LAD, week == 12 & output == val_data_LAD_reps$run$output[k]))
plot5 <- PlotReplicates(val_het_LAD_reps$samples[,1:100,k], locs = which(val_data_LAD_reps$location %in% LAD_SE), Replicates = rep_data$data) + 
  labs(x = 'LAD', y = 'Deaths')
k <- 6;rep_data <- ProcessData(subset(output_LAD, week == 12 & output == val_data_LAD_reps$run$output[k]))
plot6 <- PlotReplicates(val_het_LAD_reps$samples[,1:100,k], locs = which(val_data_LAD_reps$location %in% LAD_SE), Replicates = rep_data$data) + 
  labs(x = 'LAD', y = 'Deaths')
cowplot::plot_grid(plot1,plot2,plot3,plot4,plot5,plot6)
dev.off()

png('plots/samples_het_SW_reps.png', width = 1344, height = 960, res = 200)
k <- 1;rep_data <- ProcessData(subset(output_LAD, week == 12 & output == val_data_LAD_reps$run$output[k]))
plot1 <- PlotReplicates(val_het_LAD_reps$samples[,1:100,k], locs = which(val_data_LAD_reps$location %in% LAD_SW), Replicates = rep_data$data) + 
  labs(x = 'LAD', y = 'Deaths')
k <- 2;rep_data <- ProcessData(subset(output_LAD, week == 12 & output == val_data_LAD_reps$run$output[k]))
plot2 <- PlotReplicates(val_het_LAD_reps$samples[,1:100,k], locs = which(val_data_LAD_reps$location %in% LAD_SW), Replicates = rep_data$data) + 
  labs(x = 'LAD', y = 'Deaths')
k <- 3;rep_data <- ProcessData(subset(output_LAD, week == 12 & output == val_data_LAD_reps$run$output[k]))
plot3 <- PlotReplicates(val_het_LAD_reps$samples[,1:100,k], locs = which(val_data_LAD_reps$location %in% LAD_SW), Replicates = rep_data$data) + 
  labs(x = 'LAD', y = 'Deaths')
k <- 4;rep_data <- ProcessData(subset(output_LAD, week == 12 & output == val_data_LAD_reps$run$output[k]))
plot4 <- PlotReplicates(val_het_LAD_reps$samples[,1:100,k], locs = which(val_data_LAD_reps$location %in% LAD_SW), Replicates = rep_data$data) + 
  labs(x = 'LAD', y = 'Deaths')
k <- 5;rep_data <- ProcessData(subset(output_LAD, week == 12 & output == val_data_LAD_reps$run$output[k]))
plot5 <- PlotReplicates(val_het_LAD_reps$samples[,1:100,k], locs = which(val_data_LAD_reps$location %in% LAD_SW), Replicates = rep_data$data) + 
  labs(x = 'LAD', y = 'Deaths')
k <- 6;rep_data <- ProcessData(subset(output_LAD, week == 12 & output == val_data_LAD_reps$run$output[k]))
plot6 <- PlotReplicates(val_het_LAD_reps$samples[,1:100,k], locs = which(val_data_LAD_reps$location %in% LAD_SW), Replicates = rep_data$data) + 
  labs(x = 'LAD', y = 'Deaths')
cowplot::plot_grid(plot1,plot2,plot3,plot4,plot5,plot6)
dev.off()



#### Comparing experiments ####
results_all <- readRDS('data/samples_final_LAD/results_all.rds')
results_ward <- readRDS('data/results_ward.rds')
results_ward$Region <- as.factor(results_ward$Region)

# Only consider hetGP for plotting
tmp <- subset(results_all[,-c(3:4)], Basis == 'None')
tmp2 <- subset(results_all[,-c(3:4)], Basis == 'PLNPCA')
tmp3 <- subset(results_all[,-c(3:4)], Basis == 'SVD')
tmp$Basis <- tmp2$Basis <- tmp3$Basis <- NULL
colnames(tmp)[3:4] <- c('In95_Single', 'RMSE_Single')
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


plot1 <- ggplot(subset(tmp, Region == 'All'), aes(RMSE_Single, RMSE_Basis, col = as.factor(n))) +
  geom_point() +
  geom_abline() +
  theme_bw() +
  scale_colour_manual(values = viridis(100)[c(11,51,91)]) +
  labs(col = '\ \ n', x = 'RMSE, single emulator', y = 'RMSE, PLNPCA-LAD')

plot2 <- ggplot(subset(tmp, Region == 'All'), aes(RMSE_Basis_Ward, RMSE_Basis, col = as.factor(n))) +
  geom_point() +
  geom_abline() +
  theme_bw() +
  scale_colour_manual(values = viridis(100)[c(11,51,91)]) +
  labs(col = '\ \ n', x = 'RMSE, PLNPCA-ward', y = 'RMSE, PLNPCA-LAD')

png('plots/total_vs_basis_all_het.png', width = 1344, height = 500, res = 200)
cowplot::plot_grid(plot1, plot2, nrow = 1)
dev.off()


png('plots/total_vs_basis_region_het.png', width = 1344, height = 960, res = 200)
ggplot(subset(tmp, !(Region == 'All')), aes(RMSE_Single, RMSE_Basis, col = as.factor(n))) +
  geom_point() +
  geom_abline() +
  theme_bw() +
  facet_wrap(vars(Region)) +
  scale_colour_manual(values = viridis(100)[c(11,51,91)]) +
  labs(col = '\ \ n', x = 'RMSE, single emulator', y = 'RMSE, PLNPCA-LAD')
dev.off()

png('plots/basis_region_het.png', width = 1344, height = 960, res = 200)
ggplot(subset(tmp, !(Region == 'All')), aes(RMSE_Basis_Ward, RMSE_Basis, col = as.factor(n))) +
  geom_point() +
  geom_abline() +
  theme_bw() +
  facet_wrap(vars(Region)) +
  scale_colour_manual(values = viridis(100)[c(11,51,91)]) +
  labs(col = '\ \ n', x = 'RMSE, PLNPCA-ward', y = 'RMSE, PLNPCA-LAD')
dev.off()

# ggplot(subset(tmp, Region == 'All'), aes(RMSE_SVD, RMSE_Basis, col = as.factor(n))) +
#   geom_point() +
#   geom_abline() +
#   theme_bw() +
#   scale_colour_manual(values = viridis(100)[c(11,51,91)]) +
#   labs(col = '\ \ n', x = 'RMSE, SVD-LAD', y = 'RMSE, PLNPCA-LAD')
# 
# ggplot(subset(tmp, !(Region == 'All')), aes(RMSE_SVD, RMSE_Basis, col = as.factor(n))) +
#   geom_point() +
#   geom_abline() +
#   theme_bw() +
#   facet_wrap(vars(Region)) +
#   scale_colour_manual(values = viridis(100)[c(11,51,91)]) +
#   labs(col = '\ \ n', x = 'RMSE, SVD-LAD', y = 'RMSE, PLNPCA-LAD')












#### Plotting totals ####
#experiments <- readRDS('data/samplesBU/experiments.rds')
experiments <- readRDS('data/samples_final_LAD/results_all.rds')

# Reformat
tmp <- subset(experiments[,-c(3:4)], Basis == 'None')
tmp2 <- subset(experiments[,-c(3:4)], Basis == 'PLNPCA')
tmp$Basis <- tmp2$Basis <- NULL
colnames(tmp)[3:4] <- c('In95_Single', 'RMSE_Single')
colnames(tmp2)[3:4] <- c('In95_Basis', 'RMSE_Basis')
tmp <- tmp %>% left_join(tmp2, by = c('n', 'seed', 'Region'))


# Comparing RMSE for total, LAD_basis vs single
ggplot(subset(tmp, Region == 'All'), aes(RMSE_Single, RMSE_Basis, col = as.factor(n))) +
  geom_point() +
  geom_abline() +
  scale_colour_manual(values = viridis(100)[c(11,51,91)]) +
  theme_bw() +
  labs(col = 'n')

# Same plot, but for ward_basis vs LAD_basis
# Plot alongside above - likely shows LAD/ward very similar, then how vs single


# Comparing RMSE for all regions
# Just do for LAD
ggplot(subset(tmp, !(Region == 'All')), aes(RMSE_Single, RMSE_Basis, col = as.factor(n))) +
  geom_point() +
  geom_abline() +
  facet_wrap(vars(Region)) +
  scale_colour_manual(values = viridis(100)[c(11,51,91)]) +
  theme_bw() +
  labs(col = 'n')











colnames(experiments)[3:6] <- c('Basis', 'Total', 'Basis, hetGP', 'Total, hetGP')
plot_data <- melt(experiments[,c(1,3,4,5,6)], id.vars = 'n')
colnames(plot_data)[2:3] <- c('Emulator', 'RMSE')
ggplot(plot_data, aes(as.factor(n), RMSE, fill = Emulator)) +
  geom_boxplot()

png('plots/total_rmse.png', width = 1344, height = 960, res = 200)
ggplot(plot_data, aes(x = RMSE, col = Emulator, linetype = Emulator)) +
  stat_density(geom = 'line', position = 'identity', size = 0.65) +
  facet_wrap(vars(n), nrow = 3) +
  theme_bw() +
  scale_colour_manual(values = viridis(100)[c(1,80,1,80)]) +
  scale_linetype_manual(values = c(1,1,2,2))
dev.off()

tmp <- experiments[c(1,7:10)]
colnames(tmp)[2:5] <- c('Basis', 'Total', 'Basis, hetGP', 'Total, hetGP')
plot_data <- melt(tmp, id.vars = 'n')
colnames(plot_data)[2:3] <- c('Emulator', 'In95')
plot_data$In95[which(plot_data$n == 100)] <- plot_data$In95[which(plot_data$n == 100)] / 150
plot_data$In95[which(plot_data$n == 150)] <- plot_data$In95[which(plot_data$n == 150)] / 100
plot_data$In95[which(plot_data$n == 200)] <- plot_data$In95[which(plot_data$n == 200)] / 50
ggplot(plot_data, aes(In95, col = Emulator, linetype = Emulator)) +
  stat_density(geom = 'line', position = 'identity', size = 0.65) +
  facet_wrap(vars(n), nrow = 3) +
  geom_vline(xintercept = 0.95) +
  theme_bw() +
  scale_colour_manual(values = viridis(100)[c(1,80,1,80)]) +
  scale_linetype_manual(values = c(1,1,2,2))








ward_samples <- AggregateSamples(val_response_Het$samples)

saveRDS(ward_samples, file = 'data/ward_samples.rds')
saveRDS(lad_samples, file = 'data/lad_samples.rds')

summary(ward_samples$mean)
summary(lad_samples$mean)
summary(ward_samples$mean - lad_samples$mean)
summary(ward_samples$lower - lad_samples$lower)
summary(ward_samples$upper - lad_samples$upper)




#### Old plots ####

# plot1 <- ggplot(experiments, aes(Basis, Total, col = as.factor(n))) +
#   geom_point() +
#   geom_abline() +
#   scale_colour_manual(values = viridis(100)[c(11,51,91)]) +
#   theme_bw() +
#   labs(col = 'n') +
#   xlim(0.29,0.98) +
#   ylim(0.29,0.98)
# plot2 <- ggplot(experiments, aes(`Basis, hetGP`, `Total, hetGP`, col = as.factor(n))) +
#   geom_point() +
#   geom_abline() +
#   scale_colour_manual(values = viridis(100)[c(11,51,91)]) +
#   theme_bw() +
#   labs(col = 'n') +
#   xlim(0.29,0.98) +
#   ylim(0.29,0.98)
# png('plots/total_vs_basis.png', width = 1344, height = 600, res = 200)
# cowplot::plot_grid(plot1, plot2, nrow=1)
# dev.off()
# 
# # Plotting regions
# experiments_E12000001 <- readRDS("data/samplesBU/experiments_E12000001.rds")
# experiments_E12000002 <- readRDS("data/samplesBU/experiments_E12000002.rds")
# experiments_E12000003 <- readRDS("data/samplesBU/experiments_E12000003.rds")
# experiments_E12000004 <- readRDS("data/samplesBU/experiments_E12000004.rds")
# experiments_E12000005 <- readRDS("data/samplesBU/experiments_E12000005.rds")
# experiments_E12000006 <- readRDS("data/samplesBU/experiments_E12000006.rds")
# experiments_E12000007 <- readRDS("data/samplesBU/experiments_E12000007.rds")
# experiments_E12000008 <- readRDS("data/samplesBU/experiments_E12000008.rds")
# experiments_E12000009 <- readRDS("data/samplesBU/experiments_E12000009.rds")
# 
# experiments_all <- rbind(data.frame(experiments_E12000001, Region = 'NE'),
#                          data.frame(experiments_E12000002, Region = 'NW'),
#                          data.frame(experiments_E12000003, Region = 'YO'),
#                          data.frame(experiments_E12000004, Region = 'EM'),
#                          data.frame(experiments_E12000005, Region = 'WM'),
#                          data.frame(experiments_E12000006, Region = 'EE'),
#                          data.frame(experiments_E12000007, Region = 'LO'),
#                          data.frame(experiments_E12000008, Region = 'SE'),
#                          data.frame(experiments_E12000009, Region = 'SW'))
# 
# colnames(experiments_all)[3:6] <- c('Basis', 'Total', 'Basis, hetGP', 'Total, hetGP')
# plot_data <- melt(experiments_all[,c(1,3,4,5,6,11)], id.vars = c('n', 'Region'))
# colnames(plot_data)[3:4] <- c('Emulator', 'RMSE')
# 
# ggplot(plot_data, aes(Region, RMSE, fill = Emulator)) +
#   geom_boxplot() +
#   facet_wrap(vars(n), nrow = 2)
# 
# ggplot(plot_data, aes(as.factor(n), RMSE, fill = Emulator)) +
#   geom_boxplot() +
#   facet_wrap(vars(Region), nrow = 2)
# 
# ggplot(subset(plot_data, n == 150), aes(Region, RMSE, fill = Emulator)) +
#   geom_boxplot() +
#   theme_bw() +
#   scale_fill_manual(values = viridis(100)[c(1,80,31,100)])
# 
# ggplot(plot_data, aes(as.factor(n), RMSE, fill = Emulator)) +
#   geom_boxplot() +
#   facet_wrap(vars(Region), nrow = 3) +
#   theme_bw() +
#   scale_fill_manual(values = viridis(100)[c(1,100,1,100)])
# 
# ggplot(experiments_all, aes(Basis, Total, col = as.factor(n))) +
#   geom_point() +
#   geom_abline() +
#   scale_colour_manual(values = viridis(100)[seq(1,100,len=3)]) +
#   facet_wrap(vars(Region)) +
#   theme_bw()
# 
# ggplot(experiments_all, aes(`Basis, hetGP`, `Total, hetGP`, col = as.factor(n))) +
#   geom_point() +
#   geom_abline() +
#   scale_colour_manual(values = viridis(100)[c(11,51,91)]) +
#   facet_wrap(vars(Region)) +
#   theme_bw() +
#   labs(col = 'n')
# 
# experiments_all_total <- rbind(experiments_all,
#                                cbind(experiments, Region = 'Total'))
# experiments_all_total$Region[which(experiments_all_total$Region != 'Total')] <- 'Region'
# 
# png('plots/total_vs_basis_region.png', width = 1344, height = 550, res = 200)
# ggplot(experiments_all_total, aes(Basis, Total, col = Region)) +
#   geom_point() +
#   geom_abline() +
#   scale_colour_manual(values = viridis(100)[c(80,1)]) +
#   facet_wrap(vars(n)) +
#   theme_bw() +
#   labs(col = 'Output', y = 'Single emulator')
# dev.off()
# 
# png('plots/total_vs_basis_region_het.png', width = 1344, height = 550, res = 200)
# ggplot(experiments_all_total, aes(`Basis, hetGP`, `Total, hetGP`, col = Region)) +
#   geom_point() +
#   geom_abline() +
#   scale_colour_manual(values = viridis(100)[c(80,1)]) +
#   facet_wrap(vars(n)) +
#   theme_bw()+
#   labs(col = 'Output', y = 'Single emulator, hetGP')
# dev.off()



