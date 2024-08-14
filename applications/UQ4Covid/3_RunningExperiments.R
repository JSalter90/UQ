# Comparing different basis + emulator approaches
# Considering different sizes of training set, with/without replicates, different regions, etc.
setwd("~/Dropbox/Exeter/Projects/UQ4Covid/CountBasis_paper")
source('code/1_ProcessData.R') # also loads in 0_CountBasis.R
setwd('~/Dropbox/UQ/')
source('code/PlotFunctions.R')
source('code/Gasp.R')
setwd("~/Dropbox/Exeter/Projects/UQ4Covid/CountBasis_paper")

week12_data <- subset(output_LAD, week == 12)

set.seed(3929)
experiments <- data.frame(n = rep(c(100,150,200), each = 50),
                          seed = sample(1:10^6, 150))

# Test
# TrainEmulator('data/test', design = design_em, data = week12_data, 
#               n = experiments$n[1], seed = experiments$seed[1], 
#               basis = 'PLNPCA', reg = 'All', hetGP = TRUE, standardGP = TRUE, reps = FALSE)
# FitMulti('data/test', experiments[1:10,], design = design_em, data = week12_data, 
#               basis = FALSE, reg = 'All', hetGP = TRUE, standardGP = TRUE, reps = FALSE)

# Count basis, no reps
FitMulti(path = 'data/samples_final_LAD/PLNPCA', experiments, design_em, week12_data,
         basis = 'PLNPCA', 
         hetGP = TRUE,
         standardGP = TRUE,
         reps = FALSE)

# Some fail initially - re-run those without output
failed <- numeric(nrow(experiments))
for (i in 1:length(failed)){
  failed[i] <- length(list.files('data/samples_final_LAD/PLNPCA', pattern = paste0(experiments$seed[i]))) == 0
}
failed 

FitMulti(path = 'data/samples_final_LAD/PLNPCA', experiments[which(failed==1),], design_em, week12_data,
         basis = 'PLNPCA', 
         hetGP = TRUE,
         standardGP = TRUE,
         reps = FALSE)


# SVD basis, no reps
FitMulti(path = 'data/samples_final_LAD/SVD', experiments, design_em, week12_data,
         basis = 'SVD',
         hetGP = TRUE,
         standardGP = TRUE,
         reps = FALSE)

# Some fail initially - re-run those without output
failed <- numeric(nrow(experiments))
for (i in 1:length(failed)){
  failed[i] <- length(list.files('data/samples_final_LAD/SVD', pattern = paste0(experiments$seed[i]))) == 0
}
failed 

FitMulti(path = 'data/samples_final_LAD/SVD', experiments[which(failed==1),], design_em, week12_data,
         basis = 'SVD',
         hetGP = TRUE,
         standardGP = TRUE,
         reps = FALSE)


# Total, no reps
FitMulti(path = 'data/samples_final_LAD/Total', experiments, design_em, week12_data,
         basis = FALSE,
         hetGP = TRUE,
         standardGP = TRUE,
         reps = FALSE,
         reg = 'All')

# All regions, no reps
regs <- c('E12000001','E12000002','E12000003','E12000004','E12000005','E12000006','E12000007','E12000008','E12000009')
for (rr in regs){
  FitMulti(path = paste0('data/samples_final_LAD/', rr), experiments, design_em, week12_data,
           basis = FALSE, hetGP = TRUE, standardGP = TRUE, reps = FALSE,
           reg = rr)
}

# Ward level
# Large files, save on external drive
setwd("/Volumes/Extreme_SSD/data/CountBasis")
# Count basis, no reps
FitMulti(path = 'samples_final_ward/PLNPCA', experiments, design_em, output_ward,
         id = 'ward',
         basis = 'PLNPCA', 
         hetGP = TRUE,
         standardGP = TRUE,
         reps = FALSE)

# Some fail initially - re-run those without output
failed <- numeric(nrow(experiments))
for (i in 1:length(failed)){
  failed[i] <- length(list.files('samples_final_ward/PLNPCA', pattern = paste0(experiments$seed[i]))) == 0
}
failed 
sum(failed)

FitMulti(path = 'samples_final_ward/PLNPCA', experiments[which(failed==1),], design_em, output_ward,
         id = 'ward',
         basis = 'PLNPCA', 
         hetGP = TRUE,
         standardGP = TRUE,
         reps = FALSE)

# Some result in large values for some wards
# Change to r = 4 to see if helps
setwd("~/Dropbox/Exeter/Projects/UQ4Covid/CountBasis_paper")
results_all <-  readRDS('data/samples_final_LAD/results_all.rds')
results_ward <-  readRDS('data/results_ward.rds')
results_ward$Region <- as.factor(results_ward$Region)
tmp <- subset(results_all[,-c(3:4)], Basis == 'None')
tmp2 <- subset(results_all[,-c(3:4)], Basis == 'PLNPCA')
tmp$Basis <- tmp2$Basis <- NULL
colnames(tmp)[3:4] <- c('In95_Single', 'RMSE_Single')
colnames(tmp2)[3:4] <- c('In95_Basis', 'RMSE_Basis')
tmp <- tmp %>% left_join(tmp2, by = c('n', 'seed', 'Region'))
tmp3 <- subset(results_ward[,-c(3:4)], Basis == 'PLNPCA_ward')
colnames(tmp3)[3:4] <- c('In95_Basis_Ward', 'RMSE_Basis_Ward')
tmp <- tmp %>% left_join(tmp3, by = c('n', 'seed', 'Region'))

inds <- which(tmp$RMSE_Basis_Ward[1:150] > 1.1*tmp$RMSE_Basis[1:150])
setwd("/Volumes/Extreme_SSD/data/CountBasis")
FitMulti(path = 'samples_final_ward/PLNPCA', experiments[inds,], design_em, output_ward,
         id = 'ward',
         basis = 'PLNPCA', 
         hetGP = TRUE,
         standardGP = TRUE,
         reps = FALSE, r = 10, q = 3)

# Instead identify which locations cause issues
samp_max <- matrix(0, 8071, 150)
for (i in 1:nrow(experiments)){
  test <- list.files('samples_final_ward/PLNPCA/', pattern = paste0(experiments$seed[i]))
  if (length(test) > 0){
    All_em <- readRDS(paste0('samples_final_ward/PLNPCA/basis_PLNPCA_em_het_', experiments$seed[i], '.rds'))
    tmp <- apply(All_em$samples, 1, max)
    samp_max[,i] <- tmp
  }
}

sort(apply(samp_max, 2, max))
sort(-apply(samp_max, 1, max))

# Remove a few, see if helps
rm_ind <- order(-apply(samp_max, 1, max))[1:8]
all_wards <- unique(output_ward$ward)
tmp <- subset(output_ward, ward %in% all_wards[rm_ind] & !(duplicated(ward)))
subset(ward2019, WD19CD %in% tmp$ward)
rm_ind <- order(-apply(samp_max, 1, max))[1:8]
rm_ward <- tmp$ward

setwd("~/Dropbox/Exeter/Projects/UQ4Covid/CountBasis_paper")
results_all <-  readRDS('data/samples_final_LAD/results_all.rds')
results_ward <-  readRDS('data/results_ward.rds')
results_ward$Region <- as.factor(results_ward$Region)
tmp <- subset(results_all[,-c(3:4)], Basis == 'None')
tmp2 <- subset(results_all[,-c(3:4)], Basis == 'PLNPCA')
tmp$Basis <- tmp2$Basis <- NULL
colnames(tmp)[3:4] <- c('In95_Single', 'RMSE_Single')
colnames(tmp2)[3:4] <- c('In95_Basis', 'RMSE_Basis')
tmp <- tmp %>% left_join(tmp2, by = c('n', 'seed', 'Region'))
tmp3 <- subset(results_ward[,-c(3:4)], Basis == 'PLNPCA_ward')
colnames(tmp3)[3:4] <- c('In95_Basis_Ward', 'RMSE_Basis_Ward')
tmp <- tmp %>% left_join(tmp3, by = c('n', 'seed', 'Region'))

inds <- which(tmp$RMSE_Basis_Ward[1:150] > 1.1*tmp$RMSE_Basis[1:150])
setwd("/Volumes/Extreme_SSD/data/CountBasis")
FitMulti(path = 'samples_final_ward/PLNPCA', experiments[inds,], design_em, subset(output_ward, !(ward %in% rm_ward)),
         id = 'ward',
         basis = 'PLNPCA', 
         hetGP = TRUE,
         standardGP = TRUE,
         reps = FALSE, r = 10, q = 3)

samp_max <- matrix(0, 8071, length(inds))
for (i in 1:length(inds)){
  test <- list.files('samples_final_ward/PLNPCA/', pattern = paste0(experiments$seed[inds[i]]))
  if (length(test) > 0){
    All_em <- readRDS(paste0('samples_final_ward/PLNPCA/basis_PLNPCA_em_het_', experiments$seed[inds[i]], '.rds'))
    tmp <- apply(All_em$samples, 1, max)
    samp_max[1:length(tmp),i] <- tmp
  }
}
sort(-apply(samp_max, 1, max))[1:20] # 3 wards clearly outliers
rm_ind2 <- order(-apply(samp_max, 1, max))[1:3]
tmp <- subset(output_ward, ward %in% all_wards[-which(all_wards %in% rm_ward[1:8])][rm_ind2] & !(duplicated(ward)))
subset(ward2019, WD19CD %in% tmp$ward)
rm_ward <- c(rm_ward, tmp$ward)
subset(output_ward, ward %in% rm_ward & !(duplicated(ward)))

# For the 3 remaining issues, different wards causing it
experiments[inds,]
tt <- readRDS("/Volumes/Extreme_SSD/data/CountBasis/samples_final_ward/PLNPCA/basis_PLNPCA_em_het_45138.rds")
tt2 <- apply(tt$samples, 1, max)
all_wards[-which(all_wards %in% rm_ward)][which(tt2 > 1e7)] # "E05011091" "E05011644"
subset(ward2019, WD19CD %in% c("E05011091","E05011644"))

tt <- readRDS("/Volumes/Extreme_SSD/data/CountBasis/samples_final_ward/PLNPCA/basis_PLNPCA_em_het_232487.rds")
tt2 <- apply(tt$samples, 1, max)
all_wards[-which(all_wards %in% rm_ward)][which(tt2 > 1e7)] # "E05008330" "E05011201" "E05012762"
subset(ward2019, WD19CD %in% c("E05008330","E05011201","E05012762"))

tt <- readRDS("/Volumes/Extreme_SSD/data/CountBasis/samples_final_ward/PLNPCA/basis_PLNPCA_em_het_945202.rds")
tt2 <- apply(tt$samples, 1, max)
all_wards[-which(all_wards %in% c(rm_ward, "E05010674"))][which(tt2 > 1e7)] # "E05010061"
subset(ward2019, WD19CD %in% c("E05010061"))

TrainEmulator('samples_final_ward/PLNPCA', design = design_em, 
              data = subset(output_ward, !(ward %in% rm_ward) & !(ward %in% c("E05011091","E05011644"))), id = 'ward',
              n = experiments$n[inds[1]], seed = experiments$seed[inds[1]],
              basis = 'PLNPCA', reg = 'All', hetGP = TRUE, standardGP = TRUE, reps = FALSE, r = 10, q = 3)

TrainEmulator('samples_final_ward/PLNPCA', design = design_em, 
              data = subset(output_ward, !(ward %in% rm_ward) & !(ward %in% c("E05008330","E05011201","E05012762"))), id = 'ward',
              n = experiments$n[inds[2]], seed = experiments$seed[inds[2]],
              basis = 'PLNPCA', reg = 'All', hetGP = TRUE, standardGP = TRUE, reps = FALSE, r = 10, q = 3)

TrainEmulator('samples_final_ward/PLNPCA', design = design_em, 
              data = subset(output_ward, !(ward %in% rm_ward) & !(ward %in% c("E05010674", "E05010061"))), id = 'ward',
              n = experiments$n[inds[3]], seed = experiments$seed[inds[3]],
              basis = 'PLNPCA', reg = 'All', hetGP = TRUE, standardGP = TRUE, reps = FALSE, r = 10, q = 3)






# Repeat for reps






