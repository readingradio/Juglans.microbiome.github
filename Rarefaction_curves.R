setwd("/Users/will1809/OneDrive - purdue.edu/TCD microbiome spring 2017")

library(tidyverse)
library(readxl)

source('baxter.R')

metadata <- get_metadata()

branch <- read_tsv('Spr2017OTU.Jan2018.opti_mcc.branch.groups.rarefaction')  %>% select(-contains("lci-"), -contains("hci-")) %>% gather(-numsampled, key=sample, value=sobs) %>% mutate(sample=str_replace_all(sample, pattern="0.03-", replacement=""))  %>% drop_na()
groups  <- read_tsv('Spr2017OTU.Jan2018.opti_mcc.groups.rarefaction')  %>% select(-contains("lci-"), -contains("hci-")) %>% gather(-numsampled, key=sample, value=sobs) %>% mutate(sample=str_replace_all(sample, pattern="0.03-", replacement=""))  %>% drop_na()
leaves <- read_tsv('Spr2017OTU.Jan2018.opti_mcc.leaves.groups.rarefaction') %>% select(-contains("lci-"), -contains("hci-")) %>% gather(-numsampled, key=sample, value=sobs) %>% mutate(sample=str_replace_all(sample, pattern="0.03-", replacement=""))  %>% drop_na()

branch$sample <- gsub("-", "_", branch$sample) %>% (function (x) gsub("_L001", "", x)) %>% (function (x) gsub("_001", "", x)) %>% (function (x) gsub("_R1", "", x))
leaves$sample <- gsub("-", "_", leaves$sample) %>% (function (x) gsub("_L001", "", x)) %>% (function (x) gsub("_001", "", x)) %>% (function (x) gsub("_R1", "", x))

# note drop_na() resulted in 10 fold reduction in number of rows

# change organization of data frames so each observation is a row

alpha.branch <- read_tsv(file="Spr2017OTU.Jan2018.opti_mcc.branch.groups.ave-std.summary",
col_types=cols(group = col_character())) %>%
filter(method=='ave') %>%
select(group, sobs, shannon, invsimpson, coverage) %>%
gather(-group, key=metric, value=value)

alpha.leaves <- read_tsv(file="Spr2017OTU.Jan2018.opti_mcc.leaves.groups.ave-std.summary",
col_types=cols(group = col_character())) %>%
filter(method=='ave') %>%
select(group, sobs, shannon, invsimpson, coverage) %>%
gather(-group, key=metric, value=value)

alpha.branch$sample <- gsub("-", "_", alpha.branch$group) %>% (function (x) gsub("_L001", "", x)) %>% (function (x) gsub("_001", "", x)) %>% (function (x) gsub("_R1", "", x))
alpha.leaves$sample <- gsub("-", "_", alpha.leaves$group) %>% (function (x) gsub("_L001", "", x)) %>% (function (x) gsub("_001", "", x)) %>% (function (x) gsub("_R1", "", x))

metadata_branch <- inner_join(metadata, branch)
metadata_leaves <- inner_join(metadata, leaves)

#ggplot(metadata_branch, aes(x=numsampled, y=sobs, group=sample)) +
#geom_line()

#ggplot(groups, aes(x=numsampled, y=sobs, group=sample)) +
#geom_line()

setdiff(branch$sample, metadata_branch$sample)
setdiff(leaves$sample, metadata_leaves$sample)
setdiff(alpha.branch$sample, metadata_branch$sample)
setdiff(alpha.leaves$sample, metadata_leaves$sample)

cbind(sort(unique(alpha.leaves$sample)), sort(unique(metadata_leaves$sample)))

theme(axis.text=element_text(size=20))

# Branches
quartz(); ggplot(metadata_branch, aes(x=numsampled, y=sobs, group=sample, color=State)) + geom_line() + labs(main="Branches", x="Sequences Sampled", y="Unique OTUs") + theme(axis.text=element_text(size=14), axis.title=element_text(size=18), legend.text=element_text(size=14), legend.title=element_text(size=14))

quartz(); ggplot(metadata_branch, aes(x=numsampled, y=sobs, group=sample, color=State)) + geom_line() + labs(main="Branches", x="Sequences Sampled", y="Unique OTUs") + theme(axis.text=element_text(size=14), axis.title=element_text(size=18), legend.text=element_text(size=14), legend.title=element_text(size=14))+ coord_cartesian(xlim=c(0,20000), ylim=c(0,300))


# Leaves
quartz(); ggplot(metadata_leaves, aes(x=numsampled, y=sobs, group=sample, color=State)) + geom_line() + coord_cartesian(xlim=c(0,500), ylim=c(0,75))+ labs(main="Leaves", x="Sequences Sampled", y="Unique OTUs") + theme(axis.text=element_text(size=14), axis.title=element_text(size=18), legend.text=element_text(size=14), legend.title=element_text(size=14))
