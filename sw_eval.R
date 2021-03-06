#########################################################################
# R script: Evaluation of the performance of SWs for variable selection #
# Creator: Erik Askov Mousing                                           #
# #######################################################################

#################################
# Setup: Load data and packages #
#################################

library(mvtnorm)
library(MuMIn)
library(vegan)
library(parallel)

data.simulation <- function(sample.size=100)
{
  # sample.size: set the number of observations in the data set
  corr_xi_xj = 0 # expected correlation coefficient between independent variables
  corr_y_x1 = runif(1, -0.05, 0.05)
  corr_y_x2 = runif(1, -0.1, 0.1)
  corr_y_x3 = runif(1, -0.6, 0.6)
  corr_y_x4 = runif(1, -0.9, 0.9)
  mat = matrix(cbind(1, corr_y_x1, corr_y_x2, corr_y_x3, corr_y_x4,
                     corr_y_x1, 1, corr_xi_xj, corr_xi_xj, corr_xi_xj,
                     corr_y_x2, corr_xi_xj, 1, corr_xi_xj, corr_xi_xj,
                     corr_y_x3, corr_xi_xj, corr_xi_xj, 1, corr_xi_xj,
                     corr_y_x4, corr_xi_xj, corr_xi_xj, corr_xi_xj, 1), nrow=5)
  dat = rmvnorm(n=sample.size, mean=c(0,0,0,0,0), sigma=mat, method="svd")
  dat = as.data.frame(dat); names(dat) <- c("y", "x1", "x2", "x3", "x4")
  dat = data.frame(y=dat$y, x1=dat$x1, x2=dat$x2, x3=dat$x3, x4=dat$x4)
  return(dat)
}

dat <- data.simulation(1000)
#pairs(dat)

# Generate 1000 simulated data sets

resultlis <- list()
for(i in 1:1000)
{
  resultlis[[i]] = data.simulation(1000)
}

dat <- resultlis[[1]]

# Generate SW0 (SW baseline distribution)

cl <- makeCluster(3)
permlis <- list()

for(i in 1:100)
{
  dat <- resultlis[[i]]
  clusterExport(cl, "dat") # add data to slave processes
  clusterEvalQ(cl=cl, library(MuMIn)) # add package to slave processes
  permutation_nb <- 1000 # number of iteration in SW0 generation
  permlis[[i]] <- parSapply(cl, 1:permutation_nb, function(i)
  {
    dat$ys <- sample(dat$y)
    reg1 <- lm(ys ~ x1 + x2 + x3 + x4, data = dat, na.action = "na.fail")
    ms1 <- dredge(reg1)
    confset1 <- get.models(ms1, subset = NA)
    avgmod1 <- model.avg(confset1)
    summary(avgmod1)$importance["x1"]
  }
  )
  cat(round(100*i/1000, 2), "% \n")
  flush.console()
}

stopCluster(cl)

# Model selection using SWs

reg <- lm(y ~ x1 + x2 + x3 + x4, data = dat, na.action = "na.fail")
ms <- dredge(reg, extra = "adjR^2")
ms

# Compare
perm_weight <- permlis[[2]]

quantile(perm_weight, probs = 0.95)
summary(model.avg(get.models(ms, subset = NA)))$importance["x1"]
summary(model.avg(get.models(ms, subset = NA)))$importance["x2"]
summary(model.avg(get.models(ms, subset = NA)))$importance["x3"]
summary(model.avg(get.models(ms, subset = NA)))$importance["x4"]

#summary(reg)


varpart(dat$y, dat$x1, dat$x2, dat$x3, dat$x4)
#plot(varpart(dat$y, dat$x1, dat$x2, dat$x3, dat$x4))

anova.cca(rda(dat$y, dat$x1, cbind(dat$x2, dat$x3, dat$x4)), step = 1000)
anova.cca(rda(dat$y, dat$x2, cbind(dat$x1, dat$x3, dat$x4)), step = 1000)
anova.cca(rda(dat$y, dat$x3, cbind(dat$x1, dat$x2, dat$x4)), step = 1000)
anova.cca(rda(dat$y, dat$x4, cbind(dat$x1, dat$x2, dat$x3)), step = 1000)

