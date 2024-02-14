# Create some plots
tDataT3 <- readRDS("papers/raikoke/data/tDataT3.rds")
tDataT5 <- readRDS("papers/raikoke/data/tDataT3.rds")
tDataT7 <- readRDS("papers/raikoke/data/tDataT3.rds")
tData_regions <- readRDS("papers/raikoke/data/tDataT3.rds")

# For validation plots, see 1_Emulation.R

# Plotting T3 total vs T7 total
tmp <- data.frame(T3 = tDataT3$LogTotal + scale_output,
                  T7 = tDataT7$LogTotal + scale_output,
                  m = as.factor(design$MET))
tmp_obs <- data.frame(T3 = log(sum(all_obs[[3]]$Ash_conc_median)),
                      T7 = log(sum(all_obs[[7]]$Ash_conc_median)))
obs_sd_log <- c(max((tmp_obs[1] - log(sum(all_obs[[3]]$Percentile_10th))) / qnorm(0.9), (log(sum(all_obs[[3]]$Percentile_90th)) - tmp_obs[1])) / qnorm(0.9),
                max((tmp_obs[2] - log(sum(all_obs[[7]]$Percentile_10th))) / qnorm(0.9), (log(sum(all_obs[[7]]$Percentile_90th)) - tmp_obs[2])) / qnorm(0.9))
tmp_obs_sd <- data.frame(T3 = c(tmp_obs$T3 + qnorm(0.005)*obs_sd_log[1], tmp_obs$T3 + qnorm(0.995)*obs_sd_log[1]),
                         T7 = c(tmp_obs$T7 + qnorm(0.005)*obs_sd_log[2], tmp_obs$T7 + qnorm(0.995)*obs_sd_log[2]))

ggplot(data = tmp, aes(x = T3, y = T7, col = m)) +
  geom_point() +
  geom_vline(data = tmp_obs, aes(xintercept = T3)) +
  geom_hline(data = tmp_obs, aes(yintercept = T7)) +
  geom_vline(data = tmp_obs_sd, aes(xintercept = T3), linetype = 'dashed') +
  geom_hline(data = tmp_obs_sd, aes(yintercept = T7), linetype = 'dashed') +
  labs(x = 'Log total ash column load, T3', y = 'Log total ash column load, T7')


# Plotting N vs S, W vs E
tmp <- data.frame(Region1 = tData_regions[[1]]$LogTotal + scale_output,
                  Region2 = tData_regions[[2]]$LogTotal + scale_output,
                  m = as.factor(design$MET),
                  Region = 'North_South')

ns_line <- 48
obsR1 <- subset(all_obs[[3]], Lat >= ns_line)
obsR2 <- subset(all_obs[[3]], Lat < ns_line)
tmp_obs <- data.frame(Region1 = log(sum(obsR1$Ash_conc_median)),
                      Region2 = log(sum(obsR2$Ash_conc_median)))
obs_sd_log <- c(max((tmp_obs$Region1 - log(sum(obsR1$Percentile_10th))) / qnorm(0.9), (log(sum(obsR1$Percentile_90th)) - tmp_obs$Region1)) / qnorm(0.9),
                max((tmp_obs$Region2 - log(sum(obsR2$Percentile_10th))) / qnorm(0.9), (log(sum(obsR2$Percentile_90th)) - tmp_obs$Region2)) / qnorm(0.9))
tmp_obs_sd <- data.frame(Region1 = c(tmp_obs$Region1 + qnorm(0.005)*obs_sd_log[1], tmp_obs$Region1 + qnorm(0.995)*obs_sd_log[1]),
                         Region2 = c(tmp_obs$Region2 + qnorm(0.005)*obs_sd_log[2], tmp_obs$Region2 + qnorm(0.995)*obs_sd_log[2]))

plot1 <- ggplot(data = tmp, aes(x = Region1, y = Region2, col = m)) +
  geom_point() +
  geom_vline(data = tmp_obs, aes(xintercept = Region1)) +
  geom_hline(data = tmp_obs, aes(yintercept = Region2)) +
  geom_vline(data = tmp_obs_sd, aes(xintercept = Region1), linetype = 'dashed') +
  geom_hline(data = tmp_obs_sd, aes(yintercept = Region2), linetype = 'dashed') +
  theme(legend.position = 'none') +
  xlim(21.25,31.85) +
  ylim(23.25,31.75) +
  labs(x = 'Log total ash column load, North', y = 'Log total ash column load, South')

tmp <- data.frame(Region1 = tData_regions[[3]]$LogTotal + scale_output,
                  Region2 = tData_regions[[4]]$LogTotal + scale_output,
                  m = as.factor(design$MET),
                  Region = 'West_East')
we_line <- mean(all_obs[[3]]$Lon)
obsR3 <- subset(all_obs[[3]], Lon <= we_line)
obsR4 <- subset(all_obs[[3]], Lon > we_line)
tmp_obs <- data.frame(Region1 = log(sum(obsR3$Ash_conc_median)),
                      Region2 = log(sum(obsR4$Ash_conc_median)))
obs_sd_log <- c(max((tmp_obs$Region1 - log(sum(obsR3$Percentile_10th))) / qnorm(0.9), (log(sum(obsR3$Percentile_90th)) - tmp_obs$Region1)) / qnorm(0.9),
                max((tmp_obs$Region2 - log(sum(obsR4$Percentile_10th))) / qnorm(0.9), (log(sum(obsR4$Percentile_90th)) - tmp_obs$Region2)) / qnorm(0.9))
tmp_obs_sd <- data.frame(Region1 = c(tmp_obs$Region1 + qnorm(0.005)*obs_sd_log[1], tmp_obs$Region1 + qnorm(0.995)*obs_sd_log[1]),
                         Region2 = c(tmp_obs$Region2 + qnorm(0.005)*obs_sd_log[2], tmp_obs$Region2 + qnorm(0.995)*obs_sd_log[2]))

plot2 <- ggplot(data = tmp, aes(x = Region1, y = Region2, col = m)) +
  geom_point() +
  geom_vline(data = tmp_obs, aes(xintercept = Region1)) +
  geom_hline(data = tmp_obs, aes(yintercept = Region2)) +
  geom_vline(data = tmp_obs_sd, aes(xintercept = Region1), linetype = 'dashed') +
  geom_hline(data = tmp_obs_sd, aes(yintercept = Region2), linetype = 'dashed') +
  xlim(21.25,31.85) +
  ylim(23.25,31.75) +
  labs(x = 'Log total ash column load, West', y = 'Log total ash column load, East')

plot_grid(plot1, plot2, nrow = 1, rel_widths = c(0.95,1))

# The 4 region split
tmp <- data.frame(Region1 = tData_regions[[5]]$LogTotal + scale_output,
                  Region2 = tData_regions[[8]]$LogTotal + scale_output,
                  m = as.factor(design$MET),
                  Region = 'NW_SW')

obsR5 <- subset(all_obs[[3]], Lat >= ns_line & Lon <= we_line)
obsR6 <- subset(all_obs[[3]], Lat >= ns_line & Lon > we_line)
obsR7 <- subset(all_obs[[3]], Lat < ns_line & Lon > we_line)
obsR8 <- subset(all_obs[[3]], Lat < ns_line & Lon <= we_line)

tmp_obs <- data.frame(Region1 = log(sum(obsR5$Ash_conc_median)),
                      Region2 = log(sum(obsR8$Ash_conc_median)))
obs_sd_log <- c(max((tmp_obs$Region1 - log(sum(obsR5$Percentile_10th))) / qnorm(0.9), (log(sum(obsR5$Percentile_90th)) - tmp_obs$Region1)) / qnorm(0.9),
                max((tmp_obs$Region2 - log(sum(obsR8$Percentile_10th))) / qnorm(0.9), (log(sum(obsR8$Percentile_90th)) - tmp_obs$Region2)) / qnorm(0.9))
tmp_obs_sd <- data.frame(Region1 = c(tmp_obs$Region1 + qnorm(0.005)*obs_sd_log[1], tmp_obs$Region1 + qnorm(0.995)*obs_sd_log[1]),
                         Region2 = c(tmp_obs$Region2 + qnorm(0.005)*obs_sd_log[2], tmp_obs$Region2 + qnorm(0.995)*obs_sd_log[2]))

plot1 <- ggplot(data = tmp, aes(x = Region1, y = Region2, col = m)) +
  geom_point() +
  geom_vline(data = tmp_obs, aes(xintercept = Region1)) +
  geom_hline(data = tmp_obs, aes(yintercept = Region2)) +
  geom_vline(data = tmp_obs_sd, aes(xintercept = Region1), linetype = 'dashed') +
  geom_hline(data = tmp_obs_sd, aes(yintercept = Region2), linetype = 'dashed') +
  theme(legend.position = 'none') +
  xlim(20,31.5) +
  ylim(22.25,31) +
  labs(x = 'Log total ash column load, NW', y = 'Log total ash column load, SW')

tmp <- data.frame(Region1 = tData_regions[[6]]$LogTotal + scale_output,
                  Region2 = tData_regions[[7]]$LogTotal + scale_output,
                  m = as.factor(design$MET),
                  Region = 'NE_SE')
tmp_obs <- data.frame(Region1 = log(sum(obsR6$Ash_conc_median)),
                      Region2 = log(sum(obsR7$Ash_conc_median)))
obs_sd_log <- c(max((tmp_obs$Region1 - log(sum(obsR6$Percentile_10th))) / qnorm(0.9), (log(sum(obsR6$Percentile_90th)) - tmp_obs$Region1)) / qnorm(0.9),
                max((tmp_obs$Region2 - log(sum(obsR7$Percentile_10th))) / qnorm(0.9), (log(sum(obsR7$Percentile_90th)) - tmp_obs$Region2)) / qnorm(0.9))
tmp_obs_sd <- data.frame(Region1 = c(tmp_obs$Region1 + qnorm(0.005)*obs_sd_log[1], tmp_obs$Region1 + qnorm(0.995)*obs_sd_log[1]),
                         Region2 = c(tmp_obs$Region2 + qnorm(0.005)*obs_sd_log[2], tmp_obs$Region2 + qnorm(0.995)*obs_sd_log[2]))

plot2 <- ggplot(data = subset(tmp, Region1 > 0), aes(x = Region1, y = Region2, col = m)) +
  geom_point() +
  geom_vline(data = tmp_obs, aes(xintercept = Region1)) +
  geom_hline(data = tmp_obs, aes(yintercept = Region2)) +
  geom_vline(data = tmp_obs_sd, aes(xintercept = Region1), linetype = 'dashed') +
  geom_hline(data = tmp_obs_sd, aes(yintercept = Region2), linetype = 'dashed') +
  xlim(20,31.5) +
  ylim(22.25,31) +
  labs(x = 'Log total ash column load, NE', y = 'Log total ash column load, SE')

plot_grid(plot1, plot2, nrow = 1, rel_widths = c(0.95,1))





