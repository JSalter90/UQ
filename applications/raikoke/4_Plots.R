# Load required packages, functions
source("applications/raikoke/0_Source.R")

# Create some plots
# For validation plots, see 1_Emulation.R
tDataT1 <- readRDS("applications/raikoke/data/tDataT1.rds")
tDataT3 <- readRDS("applications/raikoke/data/tDataT3.rds")
tDataT5 <- readRDS("applications/raikoke/data/tDataT5.rds")
tData_regions <- readRDS("applications/raikoke/data/tData_regions.rds")

# Find 99% intervals
obs$Lower <- obs$Mean - qnorm(0.995)*sqrt(obs$Var)
obs$Upper <- obs$Mean + qnorm(0.995)*sqrt(obs$Var)

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
tmp_mod <- lm(T5 ~ T1, data = plot_data)

ggplot(data = subset(plot_data, m %in% 1:16), aes(x = T1, y = T5)) +
  geom_point(col = viridis(100)[20]) +
  facet_wrap(vars(m)) +
  geom_vline(data = obsT1, aes(xintercept = Mean)) +
  geom_hline(data = obsT5, aes(yintercept = Mean)) +
  geom_abline(intercept = tmp_mod$coefficients[1], slope = tmp_mod$coefficients[2], alpha = 0.5) +
  geom_vline(data = obsT1, aes(xintercept = Lower), linetype = 'dashed') +
  geom_vline(data = obsT1, aes(xintercept = Upper), linetype = 'dashed') +
  geom_hline(data = obsT5, aes(yintercept = Lower), linetype = 'dashed') +
  geom_hline(data = obsT5, aes(yintercept = Upper), linetype = 'dashed') +
  labs(x = 'Log total ash column load, T1', y = 'Log total ash column load, T5')

# Plotting N vs S, W vs E
plot_data <- data.frame(Region1 = tData_regions[[1]]$LogTotal + scale_output,
                        Region2 = tData_regions[[2]]$LogTotal + scale_output,
                        m = as.factor(design$MET),
                        Region = 'North_South')
obsR1 <- subset(obs, Type == 'N')
obsR2 <- subset(obs, Type == 'S')

tmp_mod <- lm(Region2 ~ Region1, data = plot_data)

ggplot(data = subset(plot_data, m %in% 1:16), aes(x = Region1, y = Region2, col = )) +
  facet_wrap(vars(m)) +
  geom_vline(data = obsR1, aes(xintercept = Mean)) +
  geom_hline(data = obsR2, aes(yintercept = Mean)) +
  geom_vline(data = obsR1, aes(xintercept = Lower), linetype = 'dashed') +
  geom_vline(data = obsR1, aes(xintercept = Upper), linetype = 'dashed') +
  geom_hline(data = obsR2, aes(yintercept = Lower), linetype = 'dashed') +
  geom_hline(data = obsR2, aes(yintercept = Upper), linetype = 'dashed') +
  geom_abline(intercept = tmp_mod$coefficients[1], slope = tmp_mod$coefficients[2], alpha = 0.5) +
  geom_point(alpha = 1, shape = 16, col = viridis(100)[20]) +
  xlim(21.25,31.85) +
  ylim(23.25,31.25) +
  theme(legend.position = 'none') +
  labs(x = 'Log total ash column load, North', y = 'Log total ash column load, South')

plot_data <- data.frame(Region1 = tData_regions[[3]]$LogTotal + scale_output,
                        Region2 = tData_regions[[4]]$LogTotal + scale_output,
                        m = as.factor(design$MET),
                        Region = 'West_East')
obsR3 <- subset(obs, Type == 'W')
obsR4 <- subset(obs, Type == 'E')
tmp_mod <- lm(Region2 ~ Region1, data = plot_data)
ggplot(data = subset(plot_data, m %in% 1:16), aes(x = Region1, y = Region2, col = )) +
  facet_wrap(vars(m)) +
  geom_vline(data = obsR3, aes(xintercept = Mean)) +
  geom_hline(data = obsR4, aes(yintercept = Mean)) +
  geom_vline(data = obsR3, aes(xintercept = Lower), linetype = 'dashed') +
  geom_vline(data = obsR3, aes(xintercept = Upper), linetype = 'dashed') +
  geom_hline(data = obsR4, aes(yintercept = Lower), linetype = 'dashed') +
  geom_hline(data = obsR4, aes(yintercept = Upper), linetype = 'dashed') +
  geom_abline(intercept = tmp_mod$coefficients[1], slope = tmp_mod$coefficients[2], alpha = 0.5) +
  geom_point(alpha = 1, shape = 16, col = viridis(100)[20]) +
  xlim(24.25,31.85) +
  ylim(23.5,31.2) +
  theme(legend.position = 'none') +
  labs(x = 'Log total ash column load, West', y = 'Log total ash column load, East')

# 4 regions
plot_data <- data.frame(Region1 = tData_regions[[5]]$LogTotal + scale_output,
                        Region2 = tData_regions[[8]]$LogTotal + scale_output,
                        m = as.factor(design$MET),
                        Region = 'NW_SW')
obsR5 <- subset(obs, Type == 'NW')
obsR8 <- subset(obs, Type == 'SW')
tmp_mod <- lm(Region2 ~ Region1, data = plot_data)
ggplot(data = subset(plot_data, m %in% 1:16), aes(x = Region1, y = Region2, col = )) +
  facet_wrap(vars(m)) +
  geom_vline(data = obsR5, aes(xintercept = Mean)) +
  geom_hline(data = obsR8, aes(yintercept = Mean)) +
  geom_vline(data = obsR5, aes(xintercept = Lower), linetype = 'dashed') +
  geom_vline(data = obsR5, aes(xintercept = Upper), linetype = 'dashed') +
  geom_hline(data = obsR8, aes(yintercept = Lower), linetype = 'dashed') +
  geom_hline(data = obsR8, aes(yintercept = Upper), linetype = 'dashed') +
  geom_abline(intercept = tmp_mod$coefficients[1], slope = tmp_mod$coefficients[2], alpha = 0.5) +
  geom_point(alpha = 1, shape = 16, col = viridis(100)[20]) +
  xlim(21.25,31.5) +
  ylim(22.45,31) +
  theme(legend.position = 'none') +
  labs(x = 'Log total ash column load, NW', y = 'Log total ash column load, SW')

plot_data <- data.frame(Region1 = tData_regions[[6]]$LogTotal + scale_output,
                        Region2 = tData_regions[[7]]$LogTotal + scale_output,
                        m = as.factor(design$MET),
                        Region = 'NE_SE')
obsR6 <- subset(obs, Type == 'NE')
obsR7 <- subset(obs, Type == 'SE')
tmp_mod <- lm(Region2 ~ Region1, data = plot_data)
ggplot(data = subset(plot_data, m %in% 1:16), aes(x = Region1, y = Region2, col = )) +
  facet_wrap(vars(m)) +
  geom_vline(data = obsR6, aes(xintercept = Mean)) +
  geom_hline(data = obsR7, aes(yintercept = Mean)) +
  geom_vline(data = obsR6, aes(xintercept = Lower), linetype = 'dashed') +
  geom_vline(data = obsR6, aes(xintercept = Upper), linetype = 'dashed') +
  geom_hline(data = obsR7, aes(yintercept = Lower), linetype = 'dashed') +
  geom_hline(data = obsR7, aes(yintercept = Upper), linetype = 'dashed') +
  geom_abline(intercept = tmp_mod$coefficients[1], slope = tmp_mod$coefficients[2], alpha = 0.5) +
  geom_point(alpha = 1, shape = 16, col = viridis(100)[20]) +
  xlim(20,31) +
  ylim(22.25,30.15) +
  theme(legend.position = 'none') +
  labs(x = 'Log total ash column load, NE', y = 'Log total ash column load, SE')





###¢ Plotting emulator samples
EmPost <- function(Preds, tData, obs, inds, s = 100){
  # For all m
  plot_data <- NULL
  for (i in inds){
    for (m in 0:17){
      plot_data <- rbind(plot_data, data.frame(MET = m,
                                               Run = i,
                                               Samples = t(rnorm(s, Preds$met[[m+1]]$Mean[i],
                                                                 Preds$met[[m+1]]$SD[i]))))
    }
  }
  
  # Overall
  plot_data_overall <- NULL
  for (i in inds){
    plot_data_overall <- rbind(plot_data_overall, data.frame(MET = 18, 
                                                             Run = i,
                                                             Samples = t(rnorm(s, Preds$overall$Mean[i],
                                                                               Preds$overall$SD[i]))))
  }
  
  plot_data_melt <- melt(plot_data, id.vars = c('MET', 'Run'))
  plot_data_melt$variable <- NULL
  
  plot_data_overall_melt <- melt(plot_data_overall, id.vars = c('MET', 'Run'))
  plot_data_overall_melt$variable <- NULL
  
  plot_data <- rbind(data.frame(plot_data_melt, Type = 'MET'),
                     data.frame(plot_data_overall_melt, Type = 'Overall'))
  
  plot1 <- ggplot(plot_data, aes(x = value, col = Type, linetype = as.factor(MET))) +
    stat_density(geom = 'line', position = 'identity') +
    facet_wrap(vars(Run), scales = 'free_y') +
    scale_linetype_manual(values = rep(1,19), guide = 'none') +
    scale_colour_manual(values = c(viridis(100)[c(70,1)])) +
    geom_vline(data = data.frame(Run = inds, value = tData$LogTotal[inds]), aes(xintercept = value)) +
    geom_vline(xintercept = obs, linetype = 2) +
    labs(col = 'Emulator', x = 'Log total ash column load') +
    xlim(25,32)
  
  return(plot1)
}

# Load (or create) emulator predictions for the 1000 ensemble members
EnsPredT1 <- readRDS("applications/raikoke/data/EnsPredT1.rds")
EnsPredT3 <- readRDS("applications/raikoke/data/EnsPredT3.rds")
EnsPredT5 <- readRDS("applications/raikoke/data/EnsPredT5.rds")

PostT1 <- EmPost(EnsPredT1, tDataT1, subset(obs, Type == 'T1')$Mean, inds = 1:9, s = 1000)
PostT1
PostT3 <- EmPost(EnsPredT3, tDataT3, subset(obs, Type == 'T3')$Mean, inds = 1:9, s = 1000)
PostT3
PostT5 <- EmPost(EnsPredT5, tDataT5, subset(obs, Type == 'T5')$Mean, inds = 1:9, s = 1000)
PostT5
