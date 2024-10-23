source('applications/UQ4Covid/1_ProcessData.R') # also loads in 0_CountBasis.R

# Also load in generic emulation/plotting functions from github.com/JSalter90/UQ
source('code/PlotFunctions.R')
source('code/Gasp.R')

# Comparing different basis + emulator approaches
# Considering different sizes of training set, with/without replicates, different regions, etc.
week12_data <- subset(output_LAD, week == 12)
set.seed(3929)
experiments <- data.frame(n = rep(c(100,150,200), each = 50),
                          seed = sample(1:10^6, 150))

# Count basis, no reps
# This code will loop over all combinations of (n, seed) in experiments, fit the chosen emulators, and write samples
FitMulti(path = 'data/samples_final_LAD/PLNPCA', experiments, design_em, week12_data,
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
# Large files - ensure have enough storage
# Count basis, no reps
FitMulti(path = 'data/samples_final_ward/PLNPCA', experiments, design_em, output_ward,
         id = 'ward',
         basis = 'PLNPCA', 
         hetGP = TRUE,
         standardGP = TRUE,
         reps = FALSE)

# For ward level output, particularly when smaller amount of training data, some wards with generally very low counts end up with large weights, resulting in very large sampled values
# Sometimes fixed by re-training, sometimes removing 1 or 2 of these wards fixes everything