# Load required packages, functions
source("applications/raikoke/0_Source.R")

# Create some plots
# For validation plots, see 1_Emulation.R

tDataT1 <- readRDS("applications/raikoke/data/tDataT1.rds")
tDataT3 <- readRDS("applications/raikoke/data/tDataT3.rds")
tDataT5 <- readRDS("applications/raikoke/data/tDataT5.rds")
tData_regions <- readRDS("applications/raikoke/data/tData_regions.rds")

# Find 95% intervals
obs$Lower <- obs$Mean - 1.96*sqrt(obs$Var)
obs$Upper <- obs$Mean + 1.96*sqrt(obs$Var)

# In general, scale down output by 0.95 when compare to observations
# All on log scale here
scale_output <- log(0.95)

# Plotting T1 total vs T5 total
plot_data <- data.frame(T1 = tDataT1$LogTotal + scale_output,
                        T5 = tDataT5$LogTotal + scale_output,
                        m = as.factor(design$MET),
                        logMER = design$logMER,
                        msigU = design$m_sigU)

obsT1 <- subset(obs, Type == 'T1')
obsT5 <- subset(obs, Type == 'T5')

ggplot(data = tmp, aes(x = T1, y = T5, col = m)) +
  geom_point() +
  geom_vline(xintercept = c(obsT1$Mean, obsT1$Lower, obsT1$Upper), linetype = c('solid', 'dashed', 'dashed')) +
  geom_hline(yintercept = c(obsT5$Mean, obsT5$Lower, obsT5$Upper), linetype = c('solid', 'dashed', 'dashed')) +
  labs(x = 'Log total ash column load, T1', y = 'Log total ash column load, T5')

# Alternative colour schemes
ggplot(data = plot_data, aes(x = T1, y = T5, col = m)) +
  geom_point() +
  scale_color_viridis_d(option = 'D') +
  geom_vline(xintercept = c(26.75,25.8,27.55), linetype = c('solid', 'dashed', 'dashed')) +
  geom_hline(yintercept = c(26.3,25.9,26.8), linetype = c('solid', 'dashed', 'dashed'))
  
ggplot(data = subset(plot_data, m %in% 1:16), aes(x = T1, y = T5, col = msigU)) +
  geom_point() +
  facet_wrap(vars(m)) +
  geom_abline(intercept = -3.871, slope = 1.070) +
  geom_vline(xintercept = 26.75, linetype = 'solid') +
  geom_vline(xintercept = c(25.8,27.55), linetype = 'dashed') +
  geom_hline(yintercept = 26.3, linetype = 'solid') +
  geom_hline(yintercept = c(25.9,26.8), linetype = 'dashed') +
  scale_color_viridis()

# Plotting N vs S, W vs E
plot_data <- data.frame(Region1 = tData_regions[[1]]$LogTotal + scale_output,
                        Region2 = tData_regions[[2]]$LogTotal + scale_output,
                        m = as.factor(design$MET),
                        Region = 'North_South')
obsR1 <- subset(obs, Type == 'R1')
obsR2 <- subset(obs, Type == 'R2')

plot1 <- ggplot(data = plot_data, aes(x = Region1, y = Region2, col = m)) +
  geom_point() +
  scale_color_viridis_d(option = 'D') +
  geom_vline(xintercept = c(obsR1$Mean, obsR1$Lower, obsR1$Upper), linetype = c('solid', 'dashed', 'dashed')) +
  geom_hline(yintercept = c(obsR2$Mean, obsR2$Lower, obsR2$Upper), linetype = c('solid', 'dashed', 'dashed')) +
  theme(legend.position = 'none') +
  xlim(21.25,31.85) +
  ylim(23.25,31.75) +
  labs(x = 'Log total ash column load, North', y = 'Log total ash column load, South')

plot_data2 <- data.frame(Region1 = tData_regions[[3]]$LogTotal + scale_output,
                         Region2 = tData_regions[[4]]$LogTotal + scale_output,
                         m = as.factor(design$MET),
                         Region = 'West_East')
obsR3 <- subset(obs, Type == 'R3')
obsR4 <- subset(obs, Type == 'R4')

plot2 <- ggplot(data = plot_data2, aes(x = Region1, y = Region2, col = m)) +
  geom_point() +
  scale_color_viridis_d(option = 'D') +
  geom_vline(xintercept = c(obsR3$Mean, obsR3$Lower, obsR3$Upper), linetype = c('solid', 'dashed', 'dashed')) +
  geom_hline(yintercept = c(obsR4$Mean, obsR4$Lower, obsR4$Upper), linetype = c('solid', 'dashed', 'dashed')) +
  xlim(21.25,31.85) +
  ylim(23.25,31.75) +
  labs(x = 'Log total ash column load, West', y = 'Log total ash column load, East')

plot_grid(plot1, plot2, nrow = 1, rel_widths = c(0.95,1))

# The 4 region split
plot_data <- data.frame(Region1 = tData_regions[[5]]$LogTotal + scale_output,
                        Region2 = tData_regions[[8]]$LogTotal + scale_output,
                        m = as.factor(design$MET),
                        Region = 'NW_SW')
obsR5 <- subset(obs, Type == 'R5')
obsR8 <- subset(obs, Type == 'R8')

plot1 <- ggplot(data = plot_data, aes(x = Region1, y = Region2, col = m)) +
  geom_point() +
  scale_color_viridis_d(option = 'D') +
  geom_vline(xintercept = c(obsR5$Mean, obsR5$Lower, obsR5$Upper), linetype = c('solid', 'dashed', 'dashed')) +
  geom_hline(yintercept = c(obsR8$Mean, obsR8$Lower, obsR8$Upper), linetype = c('solid', 'dashed', 'dashed')) +
  theme(legend.position = 'none') +
  xlim(20,31.5) +
  ylim(22.25,31) +
  labs(x = 'Log total ash column load, NW', y = 'Log total ash column load, SW')

plot_data2 <- data.frame(Region1 = tData_regions[[6]]$LogTotal + scale_output,
                         Region2 = tData_regions[[7]]$LogTotal + scale_output,
                         m = as.factor(design$MET),
                         Region = 'NE_SE')
obsR6 <- subset(obs, Type == 'R6')
obsR7 <- subset(obs, Type == 'R7')

plot2 <- ggplot(data = subset(plot_data2, Region1 > 0), aes(x = Region1, y = Region2, col = m)) +
  geom_point() +
  scale_color_viridis_d(option = 'D') +
  geom_vline(xintercept = c(obsR6$Mean, obsR6$Lower, obsR6$Upper), linetype = c('solid', 'dashed', 'dashed')) +
  geom_hline(yintercept = c(obsR7$Mean, obsR7$Lower, obsR7$Upper), linetype = c('solid', 'dashed', 'dashed')) +
  xlim(20,31.5) +
  ylim(22.25,31) +
  labs(x = 'Log total ash column load, NE', y = 'Log total ash column load, SE')

plot_grid(plot1, plot2, nrow = 1, rel_widths = c(0.95,1))

