setwd("~/Dropbox/Exeter/Projects/UQ4Covid/CountBasis_paper")
source('code/0_CountBasis.R')

# setwd('~/Dropbox/UQ')
# source('code/Gasp.R')
# setwd('~/Dropbox/UQ')
# source('code/PlotFunctions.R')
# setwd("~/Dropbox/Exeter/Projects/UQ4Covid/CountBasis_paper")

#### Design ####
design <- read.csv('data/design.csv')
design$lock_1_restrict <- NULL
design$lock_2_release <- NULL

# Scale inputs to [-1,1]
parRanges <- read.csv('data/parRanges.csv')
parRanges <- rbind(parRanges,
                   data.frame(parameter = c('alphaTH', 'etaTH', 'alphaEP', 'alphaI1D', 'alphaHD', 'alphaI1H', 'eta'),
                              lower = c(-1.25, 0.005, -4.35, -20, -20, -4.6, 0.005),
                              upper = c(1.65, 0.085, -1.25, -2, -1.4, -1, 0.045)))

design_em <- design
for (i in 1:15){
  ind <- which(parRanges$parameter == colnames(design[i]))
  design_em[,i] <- design_em[,i] - parRanges$lower[ind]
  design_em[,i] <- design_em[,i] / ((parRanges$upper[ind] - parRanges$lower[ind])/2)
  design_em[,i] <- design_em[,i] - 1
}


#### Outputs, LAD level ####
output_LAD <- readRDS('data/output_LAD.rds')

# Grouping LADs by region
LAD_NE <- unique(subset(output_LAD, region == 'E12000001')$LAD19CD)
LAD_NW <- unique(subset(output_LAD, region == 'E12000002')$LAD19CD)
LAD_YO <- unique(subset(output_LAD, region == 'E12000003')$LAD19CD)
LAD_EM <- unique(subset(output_LAD, region == 'E12000004')$LAD19CD)
LAD_WM <- unique(subset(output_LAD, region == 'E12000005')$LAD19CD)
LAD_EE <- unique(subset(output_LAD, region == 'E12000006')$LAD19CD)
LAD_LO <- unique(subset(output_LAD, region == 'E12000007')$LAD19CD)
LAD_SE <- unique(subset(output_LAD, region == 'E12000008')$LAD19CD)
LAD_SW <- unique(subset(output_LAD, region == 'E12000009')$LAD19CD)

#' Processing model output into a matrix
#'
#' @param data Model output (number of deaths etc., by week/LAD/region etc.)
#' @param by_id Which regional aggregation to use
#' @param output Which output to use
#' 
#' @return \item{data}{matrix of data, where rows are locations, columns are different simulations}
#' \item{location}{Names/IDs for each location}
#' \item{run}{IDs and replicate number for each simulation}
#' 
#' @export
ProcessData <- function(data, by_id = 'LAD19CD', output = 'deaths'){
  # Info on locations
  location_ids <- unique(data[,by_id])
  ell <- length(location_ids)
  
  # Info on runs
  run_ids <- unique(data[,1:2]) # need combination of run id and duplicate number
  n <- nrow(run_ids)
  
  # Is the provided dataset complete?
  stopifnot(n == nrow(data) / ell)
  
  data_matrix <- matrix(0, ell, n)
  for (i in 1:n){
    tmp_data <- subset(data, output == run_ids$output[i] & replicate == run_ids$replicate[i])[,c(by_id,output)]
    # Map to same ordering as in location_ids
    tmp_data <- left_join(as.data.frame(location_ids), 
                          tmp_data, 
                          by = c('location_ids' = paste0(by_id)))
    # Add to matrix
    data_matrix[,i] <- tmp_data[,output]
  }
  
  return(list(data = data_matrix,
              location = location_ids,
              run = run_ids))
}

# # Example
# output_LAD <- readRDS('data/output_LAD.rds')
# ens_data <- ProcessData(subset(output_LAD, week == 12 & replicate == 1))
# 
# # Some checks that worked
# tmp1 <- apply(ens_data$data, 2, sum) # total deaths per run
# tmp2 <- aggregate(deaths ~ output, subset(output_LAD, week == 12 & replicate == 1), sum)$deaths
# summary(tmp1 - tmp2)
# 
# tmp1 <- apply(ens_data$data, 1, sum) # total deaths per LAD across 250 runs
# tmp2 <- aggregate(deaths ~ LAD19CD, subset(output_LAD, week == 12 & replicate == 1), sum)$deaths
# summary(tmp1 - tmp2)
# 
# # Hospitalisations instead
# ens_data <- ProcessData(subset(output_LAD, week == 12 & replicate == 1), output = 'cumH')
# tmp1 <- apply(ens_data$data, 2, sum)
# tmp2 <- aggregate(cumH ~ output, subset(output_LAD, week == 12 & replicate == 1), sum)$cumH
# summary(tmp1 - tmp2)
# 
# tmp1 <- apply(ens_data$data, 1, sum)
# tmp2 <- aggregate(cumH ~ LAD19CD, subset(output_LAD, week == 12 & replicate == 1), sum)$cumH
# summary(tmp1 - tmp2)

# Central longitude/latitude of LADs
lad2019 <- read.csv('data/LAD2019.csv') # contains NI and Scotland, model output doesn't
LonLat <- data.frame(LAD2019CD = lad2019$lad19cd,
                     Longitude = lad2019$long,
                     Latitude = lad2019$lat)

#### Outputs, ward level ####
# Full output stored on external drive, contains all wards by week x age, large file
# Process data
# setwd("/Volumes/Extreme_SSD/data/UQ4Covid/waves")
# library(RSQLite)
# con <- dbConnect(SQLite(), 'wave04_fixedseeds_summaries_default.db')
# week12 <- dbGetQuery(con, "SELECT * FROM compact WHERE week = 12")
# dbDisconnect(con)
# dim(week12) # 45197600x9
# 
# # Aggregate by age
# week12$deaths <- week12$cumHD + week12$cumCD
# week12$Hprev_mn <- NULL
# week12$week <- NULL
# week12$replicate[is.na(week12$replicate)] <- 1
# week12$age <- as.factor(week12$age)
# week12$output <- as.factor(week12$output)
# week12_sum <- aggregate(cbind(cumH,cumHD,cumCD,deaths) ~ output + replicate + ward, data = week12, sum)
# dim(week12_sum) # 5649700 = 8071x700
# 
# # Assign various spatial regions to ward level data
# setwd("~/Dropbox/Exeter/Projects/UQ4Covid/CountBasis_paper")
# wards <- read.csv("data/Ward19_Lookup.csv")
# week12_sum$ward <- factor(wards$WD19CD[week12_sum$ward], levels = wards$WD19CD)
# 
# # Also merge in LAD as have this anyway
# week12_sum <- week12_sum %>% left_join(wards[,c('WD19CD', 'LAD19CD')], by = c('ward' = 'WD19CD'))
# 
# # Now merge in region
# load('data/ward_to_region.RData') # contains ward IDs, region IDs
# week12_sum <- week12_sum %>% left_join(ward_to_region, by = c('ward' = 'ward'))
# 
# saveRDS(week12_sum, file = 'data/output_ward.rds')
# 
# # Check that data is consistent
# output_LAD <- readRDS("data/output_LAD.rds")
# output_LAD <- subset(output_LAD, week == 12)
# week12_region <- aggregate(cbind(cumH,cumHD,cumCD,deaths) ~ output + replicate + region, data = week12_sum, FUN = sum)
# week12_lad <- aggregate(cbind(cumH,cumHD,cumCD,deaths) ~ output + replicate + LAD19CD, data = week12_sum, FUN = sum)
# 
# head(week12_lad)
# head(output_LAD)
# summary(week12_lad$deaths)
# summary(output_LAD$deaths)
# summary(week12_lad$deaths - output_LAD$deaths)
# summary(week12_lad$cumH - output_LAD$cumH)

# Load in processed ward data
output_ward <- readRDS('data/output_ward.rds')

# Grouping wards by region
ward_NE <- unique(subset(output_ward, region == 'E12000001')$ward)
ward_NW <- unique(subset(output_ward, region == 'E12000002')$ward)
ward_YO <- unique(subset(output_ward, region == 'E12000003')$ward)
ward_EM <- unique(subset(output_ward, region == 'E12000004')$ward)
ward_WM <- unique(subset(output_ward, region == 'E12000005')$ward)
ward_EE <- unique(subset(output_ward, region == 'E12000006')$ward)
ward_LO <- unique(subset(output_ward, region == 'E12000007')$ward)
ward_SE <- unique(subset(output_ward, region == 'E12000008')$ward)
ward_SW <- unique(subset(output_ward, region == 'E12000009')$ward)

# Central longitude/latitude of wards
# load('data/ward_lookup_lon_lat.RData')
# ward_lookup$lon <- ward_lookup$lon/10^5 - 5.4
# ward_lookup$lat <- ward_lookup$lat/10^5 + 49.7
# saveRDS(ward_lookup, file = 'data/ward_lookup_lon_lat.rds')
ward2019 <- readRDS('data/ward_lookup_lon_lat.rds')
LonLatWard <- data.frame(ward = ward2019$WD19CD,
                         Longitude = ward2019$lon,
                         Latitude = ward2019$lat)
