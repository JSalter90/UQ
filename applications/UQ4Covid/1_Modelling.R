source('code/0_Source.R')

# Model inputs
design <- read.csv('data/design.csv')

# Corresponding outputs, on LAD level
output_LAD <- readRDS('data/output_LAD.rds')
nrow(output_LAD) # 1661100 = 339 LADs x 7 weeks x 700 simulations

# Plot some summaries of output, e.g. total hospitalisations vs deaths in week 12
output_total <- aggregate(cbind(cumH, deaths) ~ output + replicate, data = subset(output_LAD, week == 12), sum)
ggplot(output_total, aes(log(cumH), log(deaths))) + geom_point()

# Testing basis calculation
tmp_data <- subset(output_LAD, week == 12 & replicate == 1)

# Process to a data matrix (default is by LAD, output = deaths)
ens_data <- ProcessData(tmp_data)
# This object contains an ell x n matrix
dim(ens_data$data) # 339x250
# And a list of the LAD id, run ids + replicate number that correspond to each row/column
head(ens_data$location)
head(ens_data$run)








#data_fit_st <- as.data.frame(st_data) * 7
#for (i in 1:ncol(st_data)){
#  colnames(data_fit_st)[i] <- paste0('X', i)
#}
#data_fit_st <- prepare_data(data_fit_st, numeric(nrow(st_data)))
#model_st <- PLNPCA(
#  Abundance ~ 1,
#  data  = data_fit_st, 
#  ranks = 1:10
#)
#model_st <- getBestModel(model_st, "ICL") 
#save(model_st, file = 'output/model_st.RData')


