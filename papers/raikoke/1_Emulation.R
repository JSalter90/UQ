# Load train/val inds
train_inds <- readRDS("papers/raikoke/data/train_inds.rds")
val_inds <- readRDS("papers/raikoke/data/val_inds.rds")

# Load input/output for each time point
# These objects contain the 1000 runs, where
# - each input has been scaled to [-1,1]
# - noise is a uniformly, randomly sampled term on [-1,1] to protect against overfitting mean functions
# - LogTotal is the log of the total ash for the given timepoint/region
tDataT3 <- readRDS("papers/raikoke/data/tDataT3.rds")
tDataT5 <- readRDS("papers/raikoke/data/tDataT3.rds")
tDataT7 <- readRDS("papers/raikoke/data/tDataT3.rds")

# Also for each region
# R1 = N, R2 = S, R3 = W, R4 = E, R5 = NW, R6 = NE, R7 = SE, R8 = SW
# When given in a list of 8 objects, this is the order
# Sometimes objects will have Rx in the name, where x relates to the above numbering
tData_regions <- readRDS("papers/raikoke/data/tDataT3.rds")



