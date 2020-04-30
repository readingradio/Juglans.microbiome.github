setwd("/Users/will1809/OneDrive - purdue.edu/TCD microbiome spring 2017")

library(tidyverse)
library(readxl)

source('baxter.R')

metadata <- get_metadata_soilmb()

groups  <- read_tsv('soil.mb.trim.contigs.good.unique.precluster.pick.pick.opti_mcc.groups.rarefaction')  %>% select(-contains("lci-"), -contains("hci-")) %>% gather(-numsampled, key=sample, value=sobs) %>% mutate(sample=str_replace_all(sample, pattern="0.03-", replacement=""))  %>% drop_na()

# 'sample' has to match

groups$sample[!(groups$sample == 'NTC_ITS_S96')] <- groups$sample[!(groups$sample == 'NTC_ITS_S96')] %>% sub('_S_ITS.+', '', ., perl=T)

# s[!(s == 'NTC_ITS_S96')] %>% sub('_S_ITS.+', '', ., perl=T)

# note drop_na() resulted in 10 fold reduction in number of rows

# change organization of data frames so each observation is a row

alpha.soil <- read_tsv(file="soil.mb.trim.contigs.good.unique.precluster.pick.pick.opti_mcc.groups.ave-std.summary",
col_types=cols(group = col_character())) %>%
filter(method=='ave') %>%
select(group, sobs, shannon, invsimpson, coverage) %>%
gather(-group, key=metric, value=value)

alpha.soil$sample <- alpha.soil$group %>% sub('_S_ITS.+', '', ., perl=T)

metadata_soil <- inner_join(metadata, groups)

#ggplot(metadata_branch, aes(x=numsampled, y=sobs, group=sample)) +
#geom_line()

#ggplot(groups, aes(x=numsampled, y=sobs, group=sample)) +
#geom_line()

missing.sample.groups <- setdiff(groups$sample, metadata$sample)
missing.sample.meta <- setdiff(metadata$sample, groups$sample)

setdiff(groups$sample, alpha.soil$sample)
setdiff(alpha.soil$sample, groups$sample)

setdiff(metadata$sample, alpha.soil$sample)
setdiff(alpha.soil$sample, metadata$sample)

# one of the Washington trees is missing WA132_RN_9

theme(axis.text=element_text(size=20))

# Rarefaction curve
quartz(); ggplot(metadata_soil, aes(x=numsampled, y=sobs, group=sample, color=State)) + geom_line() + labs(main="Branches", x="Sequences Sampled", y="Unique OTUs") + theme(axis.text=element_text(size=14), axis.title=element_text(size=18), legend.text=element_text(size=14), legend.title=element_text(size=14))
