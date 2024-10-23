source('applications/UQ4Covid/1_ProcessData.R') # also loads in 0_CountBasis.R

#### Aggregated validation plots ####
# PLNPCA-LAD
val_data_LAD <- readRDS("applications/UQ4Covid/data/val_data_LAD.rds") # validation data
val_het_LAD <- readRDS("applications/UQ4Covid/data/val_het_LAD.rds") # emulator samples for validation set

plot1 <- ValidateSum(val_het_LAD$samples, val_data_LAD$data) + labs(title = '')
plot2 <- ValidateSum(val_het_LAD$samples, val_data_LAD$data, 
            locs = which(val_data_LAD$location %in% LAD_NE)) + labs(title = '')

# PLNPCA-ward
val_data_ward <- readRDS("applications/UQ4Covid/data/val_data_ward.rds") # validation data
# Not on GitHub (300+Mb), needs to be created locally in 2_WorkedExamples
val_het_ward <- readRDS("applications/UQ4Covid/data/val_het_ward.rds") # emulator samples for validation set

plot3 <- ValidateSum(val_het_ward$samples, val_data_ward$data) + labs(title = '')
plot4 <- ValidateSum(val_het_ward$samples, val_data_ward$data, 
                     locs = which(val_data_ward$location %in% ward_NE)) + labs(title = '')
cowplot::plot_grid(plot1, plot3, plot2, plot4)

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
cowplot::plot_grid(plot1, plot2, plot5, plot6, rel_heights = c(1,0.85))


# Repeating for other regions
# NW, YO
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
cowplot::plot_grid(plot1, plot2, plot3, plot4, rel_heights = c(1,0.85))

# EM, WM
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
cowplot::plot_grid(plot1, plot2, plot3, plot4, rel_heights = c(1,0.85))

# EE, LO
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
cowplot::plot_grid(plot1, plot2, plot3, plot4, rel_heights = c(1,0.85))

# SE, SW
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
cowplot::plot_grid(plot1, plot2, plot3, plot4, rel_heights = c(1,0.85))


#### With reps ####
# Repeats abov validation plots, but now with replicates
# All files created in 2_WorkedExamples, not stored on GitHub when large
val_data_LAD_reps <- readRDS("applications/UQ4Covid/data/val_data_LAD_reps.rds") # validation data
val_het_LAD_reps <- readRDS("applications/UQ4Covid/data/val_het_LAD_reps.rds") # emulator samples for validation set

# Total, NE
plot1 <- ValidateSum(val_het_LAD_reps$samples, val_data_LAD_reps$data) + labs(title = '')
plot2 <- ValidateSum(val_het_LAD_reps$samples, val_data_LAD_reps$data, 
                     locs = which(val_data_LAD_reps$location %in% LAD_NE)) + labs(title = '')
cowplot::plot_grid(plot1, plot2)

# NW, YO
plot1 <- ValidateSum(val_het_LAD_reps$samples, val_data_LAD_reps$data, 
                     locs = which(val_data_LAD_reps$location %in% LAD_NW)) + labs(title = '')
plot2 <- ValidateSum(val_het_LAD_reps$samples, val_data_LAD_reps$data, 
                     locs = which(val_data_LAD_reps$location %in% LAD_YO)) + labs(title = '')
cowplot::plot_grid(plot1, plot2)

# EM, WM
plot1 <- ValidateSum(val_het_LAD_reps$samples, val_data_LAD_reps$data, 
                     locs = which(val_data_LAD_reps$location %in% LAD_EM)) + labs(title = '')
plot2 <- ValidateSum(val_het_LAD_reps$samples, val_data_LAD_reps$data, 
                     locs = which(val_data_LAD_reps$location %in% LAD_WM)) + labs(title = '')
cowplot::plot_grid(plot1, plot2)

# EE, LO
plot1 <- ValidateSum(val_het_LAD_reps$samples, val_data_LAD_reps$data, 
                     locs = which(val_data_LAD_reps$location %in% LAD_EE)) + labs(title = '')
plot2 <- ValidateSum(val_het_LAD_reps$samples, val_data_LAD_reps$data, 
                     locs = which(val_data_LAD_reps$location %in% LAD_LO)) + labs(title = '') +
  scale_x_log10(breaks = c(0.1,1,10,100,1000,10000,100000,1000000), labels = c('0','1','10','100','1000','10000','','1000000')) 
cowplot::plot_grid(plot1, plot2)

# SE, SW
plot1 <- ValidateSum(val_het_LAD_reps$samples, val_data_LAD_reps$data, 
                     locs = which(val_data_LAD_reps$location %in% LAD_SE)) + labs(title = '') +
  scale_x_log10(breaks = c(0.1,1,10,100,1000,10000,100000,1000000), labels = c('0','1','10','100','1000','10000','','1000000')) 
plot2 <- ValidateSum(val_het_LAD_reps$samples, val_data_LAD_reps$data, 
                     locs = which(val_data_LAD_reps$location %in% LAD_SW)) + labs(title = '')
cowplot::plot_grid(plot1, plot2)


#### Plotting samples ####
# LAD level, by region
PlotSamplesCount(val_het_LAD$samples, locs = which(val_data_LAD$location %in% LAD_NE), runs = 1:9, Truth = val_data_LAD$data) + 
  labs(x = 'LAD', y = 'Deaths')

PlotSamplesCount(val_het_LAD$samples, locs = which(val_data_LAD$location %in% LAD_NW), runs = 1:9, Truth = val_data_LAD$data) + 
  labs(x = 'LAD', y = 'Deaths')

PlotSamplesCount(val_het_LAD$samples, locs = which(val_data_LAD$location %in% LAD_YO), runs = 1:9, Truth = val_data_LAD$data) + 
  labs(x = 'LAD', y = 'Deaths')

PlotSamplesCount(val_het_LAD$samples, locs = which(val_data_LAD$location %in% LAD_EM), runs = 1:9, Truth = val_data_LAD$data) + 
  labs(x = 'LAD', y = 'Deaths')

PlotSamplesCount(val_het_LAD$samples, locs = which(val_data_LAD$location %in% LAD_WM), runs = 1:9, Truth = val_data_LAD$data) + 
  labs(x = 'LAD', y = 'Deaths')

PlotSamplesCount(val_het_LAD$samples, locs = which(val_data_LAD$location %in% LAD_EE), runs = 1:9, Truth = val_data_LAD$data) + 
  labs(x = 'LAD', y = 'Deaths')

PlotSamplesCount(val_het_LAD$samples, locs = which(val_data_LAD$location %in% LAD_LO), runs = 1:9, Truth = val_data_LAD$data) + 
  labs(x = 'LAD', y = 'Deaths')

PlotSamplesCount(val_het_LAD$samples, locs = which(val_data_LAD$location %in% LAD_SE), runs = 1:9, Truth = val_data_LAD$data) + 
  labs(x = 'LAD', y = 'Deaths')

PlotSamplesCount(val_het_LAD$samples, locs = which(val_data_LAD$location %in% LAD_SW), runs = 1:9, Truth = val_data_LAD$data) + 
  labs(x = 'LAD', y = 'Deaths')


# Ward-level
PlotSamplesCount(val_het_ward$samples[,1:100,], locs = which(val_data_ward$location %in% ward_NE)[1:100], runs = 1:9, Truth = val_data_ward$data) + 
  labs(x = 'Ward', y = 'Deaths')

PlotSamplesCount(val_het_ward$samples[,1:100,], locs = which(val_data_ward$location %in% ward_NW)[1:100], runs = 1:9, Truth = val_data_ward$data) + 
  labs(x = 'Ward', y = 'Deaths')

PlotSamplesCount(val_het_ward$samples[,1:100,], locs = which(val_data_ward$location %in% ward_YO)[1:100], runs = 1:9, Truth = val_data_ward$data) + 
  labs(x = 'Ward', y = 'Deaths')

PlotSamplesCount(val_het_ward$samples[,1:100,], locs = which(val_data_ward$location %in% ward_EM)[1:100], runs = 1:9, Truth = val_data_ward$data) + 
  labs(x = 'Ward', y = 'Deaths')

PlotSamplesCount(val_het_ward$samples[,1:100,], locs = which(val_data_ward$location %in% ward_WM)[1:100], runs = 1:9, Truth = val_data_ward$data) + 
  labs(x = 'Ward', y = 'Deaths')

PlotSamplesCount(val_het_ward$samples[,1:100,], locs = which(val_data_ward$location %in% ward_EE)[1:100], runs = 1:9, Truth = val_data_ward$data) + 
  labs(x = 'Ward', y = 'Deaths')

PlotSamplesCount(val_het_ward$samples[,1:100,], locs = which(val_data_ward$location %in% ward_LO)[1:100], runs = 1:9, Truth = val_data_ward$data) + 
  labs(x = 'Ward', y = 'Deaths')

PlotSamplesCount(val_het_ward$samples[,1:100,], locs = which(val_data_ward$location %in% ward_SE)[1:100], runs = 1:9, Truth = val_data_ward$data) + 
  labs(x = 'Ward', y = 'Deaths')

PlotSamplesCount(val_het_ward$samples[,1:100,], locs = which(val_data_ward$location %in% ward_SW)[1:100], runs = 1:9, Truth = val_data_ward$data) + 
  labs(x = 'Ward', y = 'Deaths')

# With replicates, can loop over regions as usual
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



#### Comparing experiments ####
results_LAD <- readRDS('applications/UQ4Covid/data/results_LAD.rds')
results_ward <- readRDS('applications/UQ4Covid/data/results_ward.rds')
results_ward$Region <- as.factor(results_ward$Region)

# Only consider hetGP for plotting
tmp <- subset(results_LAD[,-c(3:4)], Basis == 'None')
tmp2 <- subset(results_LAD[,-c(3:4)], Basis == 'PLNPCA')
tmp3 <- subset(results_LAD[,-c(3:4)], Basis == 'SVD')
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
cowplot::plot_grid(plot1, plot2, nrow = 1)

# By region
ggplot(subset(tmp, !(Region == 'All')), aes(RMSE_Single, RMSE_Basis, col = as.factor(n))) +
  geom_point() +
  geom_abline() +
  theme_bw() +
  facet_wrap(vars(Region)) +
  scale_colour_manual(values = viridis(100)[c(11,51,91)]) +
  labs(col = '\ \ n', x = 'RMSE, single emulator', y = 'RMSE, PLNPCA-LAD')

# PLNPCA-LAD vs PLNPCA-ward
ggplot(subset(tmp, !(Region == 'All')), aes(RMSE_Basis_Ward, RMSE_Basis, col = as.factor(n))) +
  geom_point() +
  geom_abline() +
  theme_bw() +
  facet_wrap(vars(Region)) +
  scale_colour_manual(values = viridis(100)[c(11,51,91)]) +
  labs(col = '\ \ n', x = 'RMSE, PLNPCA-ward', y = 'RMSE, PLNPCA-LAD')



#### Plotting correlations ####
basis_deaths <- readRDS('applications/UQ4Covid/data/basis_deaths_example.rds')
val_data_LAD <- readRDS("applications/UQ4Covid/data/val_data_LAD.rds") # validation data
val_het_LAD <- readRDS("applications/UQ4Covid/data/val_het_LAD.rds") # emulator samples for validation set

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
  ctrain <- csample <- numeric(339)
  for (i in 1:339){
    ctrain[i] <- cor(basis_deaths$Data[k[j],], basis_deaths$Data[i,])
    csample[i] <- cor(c(val_het_LAD$samples[k[j],,]), c(val_het_LAD$samples[i,,]))
  }
  inds <- order(-ctrain)
  plot_data <- rbind(plot_data, data.frame(k = regions[j],
                                           LAD = 1:339,
                                           Training = ctrain[inds],
                                           Emulator = csample[inds]))
}
plot_data$k <- factor(plot_data$k, levels = regions)
plot_data <- melt(plot_data, id.vars = c('LAD', 'k'))
ggplot(plot_data, aes(LAD, value, col = variable)) + 
  geom_line(size = 0.5) +
  geom_line(data = subset(plot_data, variable == 'Training'), size = 0.75) + 
  facet_wrap(vars(k)) +
  scale_colour_manual(values = viridis(100)[c(71,21)]) +
  theme_bw() +
  labs(y = 'Correlation', col = 'Type')
