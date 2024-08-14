# How correlations pass through
# Outputs are unrelated as var[c(x)] -> 0, e.g. if 
cor(rpois(1000,1), rpois(1000,1))
cor(rpois(1000,10), rpois(1000,1000))

# But this is only true for later vectors (i.e., those that represent random noise)
# If this were true for leading vectors, suggests model is not stochastic

# If 2 outputs are perfectly correlated on basis, maintain a high correlation
cor(rpois(1000,1:100), rpois(1000,1:100))
plot(rpois(1000,1:100), rpois(1000,1:100))

# Similarly if negatively correlated
cor(rpois(1000,1:100), rpois(1000,100:1))
plot(rpois(1000,1:100), rpois(1000,100:1))

# Slightly stronger correlation maintained on a log scale
cor(log(rpois(1000,1:100)+1), log(rpois(1000,1:100)+1))

# Less correlation maintained as variance decreases
tt <- rnorm(1000,500,200)
cor(1000 + tt, 2000 + tt)
cor(rpois(1000, 1000 + tt), rpois(1000, 2000 + tt))

tt <- rnorm(1000,500,100)
cor(1000 + tt, 2000 + tt)
cor(rpois(1000, 1000 + tt), rpois(1000, 2000 + tt))

tt <- rnorm(1000,500,10)
cor(1000 + tt, 2000 + tt)
cor(rpois(1000, 1000 + tt), rpois(1000, 2000 + tt))

# But the variance is a tuning parameter
# To ensure get correlation on Poisson level, increases variance on latent level
# If low variance patterns were in fact important for correlation, would have high variance


# What happens if train a model where 2 outputs are perfectly correlated?
train_labels <- design$output[1:50]
train_data <- subset(output_LAD, week == 12 & replicate == 1 & output %in% train_labels) # 67800 = 200 runs x 339 LADs
train_data <- ProcessData(train_data)
dim(train_data$data)
cor(log(train_data$data[1,]+1), log(train_data$data[2,]+1))
tmp_basis <- CountBasis(train_data$data, rank = 10)

# Make 2 outputs identical
train_data$data[2,] <- train_data$data[1,]
tmp_basis2 <- CountBasis(train_data$data, rank = 10)

tmp_basis$tBasis[1:2,]
tmp_basis2$tBasis[1:2,]

tt1 <- tmp_basis$EnsembleMean + tmp_basis$LatentMean + tmp_basis$tBasis %*% t(tmp_basis$Coeffs)
tt2 <- tmp_basis2$EnsembleMean + tmp_basis2$LatentMean + tmp_basis2$tBasis %*% t(tmp_basis2$Coeffs)

tt1[1:2,]
tt2[1:2,]

# If sample across the ensemble - very correlated (because have high variance)
cor(rpois(1000, exp(tt1[1,])), rpois(1000, exp(tt1[2,])))
cor(rpois(1000, exp(tt2[1,])), rpois(1000, exp(tt2[2,])))

plot(rpois(1000, exp(tt1[1,])), rpois(1000, exp(tt1[2,])))

# If sample for particular run - not correlated (because assumiung zero variance)
cor(rpois(1000, exp(tt1[1,1])), rpois(1000, exp(tt1[2,1])))
cor(rpois(1000, exp(tt2[1,1])), rpois(1000, exp(tt2[2,1])))

plot(rpois(1000, exp(tt1[1,1])), rpois(1000, exp(tt1[2,1])))


cor(rpois(1000, 5 + rnorm(1000,0.01,1)), rpois(1000, 10 + rnorm(1000,0.01,1)))



# What if there's no structured variability in the data?
train_data <- subset(output_LAD, week == 12 & replicate == 1 & output %in% train_labels) # 67800 = 200 runs x 339 LADs
train_data3 <- ProcessData(train_data)
for (j in 2:50){
  train_data3$data[,j] <- train_data3$data[,1] + rpois(339,1)
}
tmp_basis3 <- CountBasis(train_data3$data, rank = 10)

apply(tmp_basis$Coeffs, 2, var)
apply(tmp_basis2$Coeffs, 2, var)
apply(tmp_basis3$Coeffs, 2, var)

tt3 <- tmp_basis3$EnsembleMean + tmp_basis3$LatentMean + tmp_basis3$tBasis %*% t(tmp_basis3$Coeffs)

# No/little correlation in the training data
cor(train_data3$data[1,], train_data3$data[2,])
cor(train_data3$data[10,], train_data3$data[200,])
cor(train_data3$data[50,], train_data3$data[52,])

# Much higher correlation in the latent structure
cor(tt3[1,], tt3[2,])
cor(tt3[10,], tt3[200,])
cor(tt3[50,], tt3[52,])

# Low correlation when sample
cor(rpois(1000, exp(tt3[1,])), rpois(1000, exp(tt3[2,])))
cor(rpois(1000, exp(tt3[10,])), rpois(1000, exp(tt3[200,])))
cor(rpois(1000, exp(tt3[50,])), rpois(1000, exp(tt3[52,])))


#### High correlation but low variability? (i.e most data identical?) ####
train_data <- subset(output_LAD, week == 12 & replicate == 1 & output %in% train_labels) # 67800 = 200 runs x 339 LADs
train_data3 <- ProcessData(train_data)
for (j in 2:50){
  train_data3$data[,j] <- train_data3$data[,1] + 100*(1:339)
}
tmp_basis3 <- CountBasis(train_data3$data, rank = 10)

apply(tmp_basis$Coeffs, 2, var)
apply(tmp_basis2$Coeffs, 2, var)
apply(tmp_basis3$Coeffs, 2, var)

tt3 <- tmp_basis3$EnsembleMean + tmp_basis3$LatentMean + tmp_basis3$tBasis %*% t(tmp_basis3$Coeffs)

cor(train_data3$data[1,], train_data3$data[2,])
cor(train_data3$data[10,], train_data3$data[200,])
cor(train_data3$data[50,], train_data3$data[52,])

cor(tt3[1,], tt3[2,])
cor(tt3[10,], tt3[200,])
cor(tt3[50,], tt3[52,])

# Low correlation when sample
cor(rpois(1000, exp(tt3[1,])), rpois(1000, exp(tt3[2,])))
cor(rpois(1000, exp(tt3[10,])), rpois(1000, exp(tt3[200,])))
cor(rpois(1000, exp(tt3[50,])), rpois(1000, exp(tt3[52,])))




#### Comparing empirical moments to theoretical ####
basis_deaths <- readRDS('data/basis_deaths_example.rds')
basis <- basis_deaths$tBasis
k <- 10
z_mean <- as.numeric(basis_deaths$EnsembleMean + basis_deaths$LatentMean + basis %*% basis_deaths$Coeffs[k,])
z_cov <- basis %*% diag(10:1) %*% t(basis)
samps <- rmvnorm(10000, z_mean, z_cov) # 1000x339
samps_p <- matrix(rpois(10000*339, exp(samps)), 10000, 339)

samps[1:10,1:6];samps_p[1:10,1:6]

# Theoretical moments
mm <- as.numeric(exp(z_mean + diag(z_cov)/2))
vv <- mm + mm^2 * (exp(diag(z_cov)) - 1)
cc <- mm[1]*mm*(exp(z_cov[1,])-1)
cc2 <- sqrt(mm[1])*sqrt(mm) + mm[1]*mm*(exp(z_cov[1,])-1)

# Sample moments
mms <- apply(samps_p, 2, mean)
plot(mm, mms)
cor(mm, mms)

vvs <- apply(samps_p, 2, var)
plot(vv, vvs)
cor(vv, vvs)

ccs <- c(cov(samps_p[,1], samps_p))
plot(cc, ccs); abline(0,1)
cor(cc, ccs)

plot(cc[-1], ccs[-1]); abline(0,1)
cor(cc[-1], ccs[-1])

plot(cc2, ccs); abline(0,1)
cor(cc2, ccs)

cor(samps_p[,1], samps_p[,2]) # sample correlation
cc[2] / sqrt(vv[1]*vv[2]) # theoretical correlation



tt <- rnorm(100,1000,0.01); tt2 <- tt+0.1
cor(tt,tt2)
cor(rpois(100,tt), rpois(100,tt2))
plot(rpois(100,tt), rpois(100,tt2))

tt <- rnorm(100,1000,100); tt2 <- tt+0.1
cor(tt,tt2)
cor(rpois(100,tt), rpois(100,tt2))
plot(rpois(100,tt), rpois(100,tt2))

cor(rpois(100,500+1:100), rpois(100,5+1:100))
plot(rpois(100,500+1:100), rpois(100,5+1:100))