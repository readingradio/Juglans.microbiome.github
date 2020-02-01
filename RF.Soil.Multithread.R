# This script performs a random forest classification on merged bacterial and fungal
# mothur data. This is in a seperate script so it can be called from shell for HPC.

# Memory and processor time issues were encountered with the newest version of
# randomForest. For this reason, the prunescale cutoff was set higher for soil
# to prevent a stack overflow. To deal with processing time issues, I used the
# packages "doParallel" and "foreach" to implement parallel processing of the
# bootstrap.

setwd("/scratch/brown/will1809/RF")
load("RF.soil.input.RData")

library(doParallel)
library(foreach)
library(tidyverse)
library(randomForest)

set.seed(579383)

registerDoParallel(24)

# %dopar% implemented with foreach implements each loop call as an independent function call
# and then parallelizes all the function calls across the number of processers set above with
# registerDoParallel(). First, we need to set up a dataframe to store each bootstrap RF result.
boot.imp.s <- data.frame(predictors=NULL, Try=NULL, MeanDecreaseGini=NULL)
boot.imp.s <- foreach (Try=1:100, .combine=rbind) %dopar% {

  soils.classify <- randomForest(response.soils ~., data = rf.data.soils, ntree = 500, proximity=T, mtry=200)
  print(soils.classify)

  # Make a data frame with predictor names and their importance
  imp.s <- importance(soils.classify)
  imp.s <- data.frame(predictors = rownames(imp.s), imp.s)

  cbind(imp.s, data.frame(Try=Try)) #)
}

# Same (and independent) as above but this time we are storing results from the
# out-of-bag error rates, so we car report the correct-incorrect classification
# rate for each group, and what types of misclassifications were made.
boot.oob <- NULL
boot.oob <- foreach (Try=1:100, .combine=rbind) %dopar% {

  soils.classify <- randomForest(response.soils ~., data = rf.data.soils, ntree = 500, proximity=T, mtry=200)
  print(soils.classify)

  # Make a data frame with predictor names and their importance
  imp.s <- importance(soils.classify)
  imp.s <- data.frame(predictors = rownames(imp.s), imp.s)

  c(boot.oob, soils.classify$err.rate[500,1])
}

save.image('RF.soil.output.RData')
quit()
