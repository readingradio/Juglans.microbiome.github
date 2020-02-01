# This script performs a random forest classification on merged bacterial and fungal
# mothur data. It plots the results as a heatmap and barchart containing bootstrapped
# importance values for predicting where a sample came from based on its microbiome
# profile with a susbset of the most important OTUs as predictor variables. It also
# outputs a taxonomy table with the identity of those OTUs.

# Memory and processor time issues were encountered with the newest version of
# randomForest. For this reason, the prunescale cutoff was set higher for soil
# to prevent a stack overflow. To deal with processing time issues, I used the
# packages "doParallel" and "foreach" to implement parallel processing of the
# bootstrap.

setwd(".")

load("R_Environments/Jnigra.microbiome.merged.WilliamsOnufrak.RData")

library(doParallel)
library(foreach)
library(tidyverse)
library(ggplot2)
library(vegan)
library(reshape2)
library(phyloseq)
library(randomForest)
library(doBy)
library(plotrix)
library(ggpubr)

#######################
##### CAULOSPHERE #####
#######################

# Set prunescale 
# Prunescale is the minimum average number of reads across samples
# that will be retained. All OTUs that do not average > prunescale
# across samples will be dropped.
prunescale.shavings = 0.05

# Prune out rare OTUs by mean relative abundance set by prunescale
tax.mean.shavings <- taxa_sums(shavings.merge)/nsamples(shavings.merge) # average number of reads for each otu
sites.prune.shavings <- prune_taxa(tax.mean.shavings > prunescale.shavings, shavings.merge)

# Make a dataframe of training data with OTUs as column and samples as rows
predictors.shavings <- t(otu_table(sites.prune.shavings))
dim(predictors.shavings)

# Make one column for our outcome/response variable 
response.shavings <- as.factor(sample_data(sites.prune.shavings)$State)

# Combine them into 1 data frame
rf.data.shavings <- data.frame(response.shavings, predictors.shavings)

# Set a random seed
set.seed(579383)

# My computer has 8 processor cores - adjust parameter of the following command as needed/desired
registerDoParallel(8)

# %dopar% implemented with foreach implements each loop call as an independent function call
# and then parallelizes all the function calls across the number of processers set above with
# registerDoParallel(). First, we need to set up a dataframe to store each bootstrap RF result.
boot.imp.s <- data.frame(predictors=NULL, Try=NULL, MeanDecreaseGini=NULL)
boot.imp.s <- foreach (Try=1:100, .combine=rbind) %dopar% {

  shavings.classify <- randomForest(response.shavings ~., data = rf.data.shavings, ntree = 500, proximity=T, mtry=200)
  print(shavings.classify)

  # Make a data frame with predictor names and their importance
  imp.s <- importance(shavings.classify)
  imp.s <- data.frame(predictors = rownames(imp.s), imp.s)

  cbind(imp.s, data.frame(Try=Try)) #)

}

# Same (and independent) as above but this time we are storing results from the
# out-of-bag error rates, so we car report the correct-incorrect classification
# rate for each group, and what types of misclassifications were made.
boot.oob <- NULL
boot.oob <- foreach (Try=1:100, .combine=rbind) %dopar% {

  shavings.classify <- randomForest(response.shavings ~., data = rf.data.shavings, ntree = 500, proximity=T, mtry=200)
  print(shavings.classify)

  # Make a data frame with predictor names and their importance
  imp.s <- importance(shavings.classify)
  imp.s <- data.frame(predictors = rownames(imp.s), imp.s)

  c(boot.oob, shavings.classify$err.rate[500,1])
}

# Look at the results
mean(boot.oob)
sd(boot.oob)
range(boot.oob)

# Create a dataframe with summary statistics from our bootstrap runs
summary.boot.s <- summaryBy(MeanDecreaseGini ~ predictors, data= boot.imp.s, FUN = c(mean, sd, std.error))

# Look at it
head(summary.boot.s)

# Order the predictor levels by importance
imp.s.sort <- arrange(summary.boot.s, desc(MeanDecreaseGini.mean))
imp.s.sort$predictors <- factor(imp.s.sort$predictors, levels = imp.s.sort$predictors)

# Select the top 40 predictors
imp.s.40 <- imp.s.sort[1:40, ]

# What are those OTUs?
# otunames.s = otunames of predictors
# r.s = A logical containing which rows from phyloseq object were recovered as top 40 predictors
otunames.s <- imp.s.40$predictors
r.s <- rownames(tax_table(shavings.merge)) %in% otunames.s

# Use r.s to pull out the taxonomic information for those OTUS in otunames.s
t.table <- as.data.frame(tax_table(shavings.merge)[r.s, ])

write.csv(t.table, "Results/RF.shavings.t.table.csv", quote=F,row.names=T)

# Reformat the taxonomic names and paste them to OTU names for axis labels
# (highest determined taxonomic rate, get rid of special characters, etc.)
axislabels <- with(t.table, paste(rownames(t.table), gsub("[sgfocp]__","",Genus,perl=T) %>% gsub("_"," ",.,perl=T) %>% gsub("[\\s]?(unclassified|ge)[\\s]?","",.,perl=T)))
names(axislabels) <- rownames(t.table)
i <- grep("(uncultured)|([01234567890-]$)+", axislabels, perl=T)
axislabels[i]<-paste(names(axislabels)[i], t.table$Family[i])

# Now we sum the abundance of each OTU by state so we can order the
# RF predictor OTUs by which state they were most abundant in, followed
# by their importance value. This makes the heatmap easier to interpret.
sums<-summaryBy(. ~ response.shavings, data = rf.data.shavings[, c('response.shavings', as.character(otunames.s))], FUN = sum, keep.names = T)
state.organize <- sapply(as.character(otunames.s), FUN = function(x) which.max(sums[,x]))

# hm is a dataframe of the abundance of the predictors, standardized and log transformed
hm <- rf.data.shavings[,as.character(otunames.s)] %>% decostand(method='max', MARGIN=2) %>% decostand(method='log')
hm$sample <- gsub('_',' ',sample_data(shavings.merge)[rownames(hm),'New.names']$New.names)

# Use the melt() function to turns the ordered dataframe into
# a new data frame that represents the value in each column as
# a new row so that each row represents an individual observation.
# This creates a dataframe with columnes for tree, the otu, and its
# standardized transformed abundance in that tree for the heatmap.
hm.melt <- melt(hm[c(order(state.organize), 41)])
names(hm.melt) <- c("Tree","Taxon","StandardizedAbundance")

# Turn predictors into a factor (re-ordered the same as the hm.melt object so they match up)
imp.s.40$predictors <- factor(imp.s.40$predictors, levels(imp.s.40$predictors)[order(state.organize)])

# Creat the heatmap with the taxon names and OTUs
plot1 <- ggplot(hm.melt, aes(x=Tree, y=Taxon, fill=StandardizedAbundance)) + geom_tile() + scale_fill_gradient(low="white", high="black") + theme(axis.text.x = element_text(angle=90, size=6), legend.position = "none", axis.text.y = element_text(hjust=0))  + scale_y_discrete(labels=axislabels)
g1 <- ggplotGrob(plot1)

# Creat the bargraph
plot2 <- ggplot(imp.s.40, aes(x = predictors, y = MeanDecreaseGini.mean)) +
  geom_bar(stat = "identity", fill = "grey") +
  ylab("Bootstrap Mean \u00b1 SD Gini Index") +
  theme_bw() +
  geom_errorbar( aes(x=predictors, ymin=MeanDecreaseGini.mean-MeanDecreaseGini.sd, ymax=MeanDecreaseGini.mean+MeanDecreaseGini.sd), width=0.75, alpha=0.9, size=0.5) +
  coord_flip() + theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
		axis.line.x=element_line(),
		axis.title=element_text(size=12),
			plot.background = element_rect(fill = NULL,colour = NA),
			panel.border = element_blank(),
			panel.grid = element_blank())
g2 <- ggplotGrob(plot2)

# Put them togethers
ggarrange(g1, g2, align='h', nrow=1, ncol=2, widths=c(4,2))

pdf(paste("Figures/Caulosphere.RF", Sys.time(), ".pdf"), 9,6)
ggarrange(g1, g2, align='h', nrow=1, ncol=2, widths=c(4,2))
dev.off()

#####################
##### BULK SOIL #####
#####################

# Set prunescale higher for big data frame
prunescale.soils = 0.2

# Prune out rare OTUs by mean relative abundance set by prunescale
tax.mean.soils <- taxa_sums(soils.merge)/nsamples(soils.merge)
sites.prune.soils <- prune_taxa(tax.mean.soils > prunescale.soils, soils.merge)
# Make a dataframe of training data with OTUs as column and samples as rows
predictors.soils <- t(otu_table(sites.prune.soils))
dim(predictors.soils)

# Make one column for our outcome/response variable 
response.soils <- as.factor(sample_data(sites.prune.soils)$State)

# Combine them into 1 data frame
rf.data.soils <- data.frame(response.soils, predictors.soils)

# To export for multithreading on HPC
save.image("R_Environments/RF.soil.input.RData")

#### IMPORTANT IMPORTANT IMPORTANT ####
#### Now run __ "RF.Soil.Multithread.R"
#### either on HPC with "RF.Soil.Multithread.R.sh"
#### with "RF.soil.input.RData" in the current folder
#### or locally with source("RF.Soil.Multithread.R")
#### It will input the data from "RF.soil.input.RData"
#### It will output the data to "RF.soil.output.RData"
#### Put the output file in the folder "R_Environments"

load("R_Environments/RF.soil.output.RData")

mean(boot.oob)
sd(boot.oob)
range(boot.oob)

summary.boot.s <- summaryBy(MeanDecreaseGini ~ predictors, data= boot.imp.s, FUN = c(mean, sd, std.error))

head(summary.boot.s)

# Order the predictor levels by importance
imp.s.sort <- arrange(summary.boot.s, desc(MeanDecreaseGini.mean))
imp.s.sort$predictors <- factor(imp.s.sort$predictors, levels = imp.s.sort$predictors)

# Select the top 10 predictors
imp.s.40 <- imp.s.sort[1:40, ]

# What are those OTUs?

otunames.s <- imp.s.40$predictors
r.s <- rownames(tax_table(soils.merge)) %in% otunames.s
t.table <- as.data.frame(tax_table(soils.merge)[r.s, ])

write.csv(t.table, "Results/RF.soil.t.table.csv", quote=F,row.names=T)

axislabels <- with(t.table, paste(rownames(t.table), gsub("[sgfocp]__","",Genus,perl=T) %>% gsub("_"," ",.,perl=T) %>% gsub("[\\s]?(unclassified|ge)[\\s]?","",.,perl=T)))
names(axislabels) <- rownames(t.table)

i <- grep("(uncultured$)|([01234567890-]+$)|(group$)", axislabels, perl=T)
axislabels[i]<-paste(names(axislabels)[i], t.table$Family[i])
j <- grep("(uncultured$)|([01234567890-]+$)|(group$)|([01234567890-]+_fa$)", axislabels, perl=T)
axislabels[j]<-paste(names(axislabels)[j], t.table$Order[j])
k <- grep("(uncultured_fa$)|([01234567890-]+_or$)", axislabels, perl=T)
axislabels[k]<-paste(names(axislabels)[k], t.table$Phylum[k])
l <- grep("[01234567890-]+$", axislabels, perl=T)
axislabels[l]<-paste(names(axislabels)[l], t.table$Class[l])

hm <- rf.data.soils[,as.character(otunames.s)] %>% decostand(method='max', MARGIN=2) %>% decostand(method='log')
hm$sample <- gsub('_',' ',sample_data(soils.merge)[rownames(hm),'New.names']$New.names)

sums<-summaryBy(. ~ response.soils, data = rf.data.soils[, c('response.soils', as.character(otunames.s))], FUN = sum, keep.names = T)
state.organize <- sapply(as.character(otunames.s), FUN = function(x) which.max(sums[,x]))

hm.melt <- melt(hm[c(order(state.organize), 41)])
names(hm.melt) <- c("Tree","Taxon","StandardizedAbundance")

imp.s.40$predictors <- factor(imp.s.40$predictors, levels(imp.s.40$predictors)[order(state.organize)])

plot1 <- ggplot(hm.melt, aes(x=Tree, y=Taxon, fill=StandardizedAbundance)) + geom_tile() + scale_fill_gradient(low="white", high="black") + theme(axis.text.x = element_text(angle=90, size=6), legend.position = "none", axis.text.y = element_text(hjust=0))  + scale_y_discrete(labels=axislabels)

#quartz()
g1 <- ggplotGrob(plot1)

plot2 <- ggplot(imp.s.40, aes(x = predictors, y = MeanDecreaseGini.mean)) +
  geom_bar(stat = "identity", fill = "grey") +
  ylab("Bootstrap Mean \u00b1 SD Gini Index") +
  theme_bw() +
  geom_errorbar( aes(x=predictors, ymin=MeanDecreaseGini.mean-MeanDecreaseGini.sd, ymax=MeanDecreaseGini.mean+MeanDecreaseGini.sd), width=0.75, alpha=0.9, size=0.5) +
  coord_flip() + theme(axis.title.y=element_blank(),
                       axis.text.y=element_blank(),
                       axis.ticks.y=element_blank(),
                       axis.line.x=element_line(),
                       axis.title=element_text(size=12),
                       plot.background = element_rect(fill = NULL,colour = NA),
                       panel.border = element_blank(),
                       panel.grid = element_blank())
g2 <- ggplotGrob(plot2)

ggarrange(g1, g2, align='h', nrow=1, ncol=2, widths=c(4,2))

pdf(paste("Figures/Soil.RF", Sys.time(), ".pdf"), 9,6)
ggarrange(g1, g2, align='h', nrow=1, ncol=2, widths=c(4,2))
dev.off()
