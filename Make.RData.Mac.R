# This script reads in everything you need to complete downstream analysis of
# the merged bacterial and fungal dataframes for soil and bacteria from J. nigra.

# It uses functions from 'baxter.nov2019.R' to read in and reformat the Excell
# metadata file(s) and then the rarefied bacterial and fungal OTU tables.

# Then it creates a phyloseq object containing all the data from the rarefied
# dataframes including samples retained in both bacterial and fungal rarefaction,
# and saves an R image for easy loading the data as a phyloseq object in the R
# environment as an .RData file

library(tidyverse)
library(phyloseq)

# Put your path to project diretory here
setwd(".")

#baxter.R has the get_metadata() function I wrote in it which is nice
source('baxter.nov2019.R')

# Put your path to mothur output diretory here
setwd("Mothur_output")

# Read in fungal shaving data and load it as a phyloseq object
fshavings.names <- (read_tsv("Phyllo.ITS2.jn.Oct2019.opti_mcc.branch.0.02.subsample.shared", col_types=cols(Group=col_character())) )$Group %>% sub('_S_ITS.+', '', ., perl=T)
fshavings.ps <- import_mothur(mothur_shared_file = "Phyllo.ITS2.jn.Oct2019.opti_mcc.branch.0.02.subsample.shared", mothur_constaxonomy_file = "phyllo.mb.its2.good.unique.precluster.pick.pick.opti_mcc.0.02.cons.taxonomy")

# number of NumOTUS to add to bark shared file column
dim(read_tsv("16s_bark_otutable_rarefied_nosingleton_fa19_ajo.shared", col_types=cols(Group=col_character()) ))[2]

# number of NumOTUS to add to soil shared file column
dim(read_tsv("16s_soil_otutable_rarefied_nosingleton_fa19_ajo.shared", col_types=cols(Group=col_character()) ))[2]

# To reformat bacterial shared and taxonomy files prior to merging with fungus,
# I did the equivalent of what follows below in BASH with 'sed'
#row.names(otu_table(bshavings.ps)) <- row.names(otu_table(bshavings.ps)) %>% sub("Otu", "BOtu", .)
#row.names(tax_table(bshavings.ps)) <- row.names(tax_table(bshavings.ps)) %>% sub("Otu", "BOtu", .)

# Also had to reformat the file to be compatible with import_mothur()
# Which consisted of:
### removing "
### adding columns name Group
### adding columns label (=0.03) and numBotus (=1775 for bark and =21433 for soil)
# First check number of OTUs
#read.csv("16s_bark_otutable_rarefied_nosingleton_fa19_ajo.shared", sep='\t') %>% dim() # 1775 OTUs
#read.csv("16s_bark_otutable_rarefied_nosingleton_fa19_ajo.shared.botu", sep='\t') %>% dim() # 1775 OTUs (1776 -1) I added a column
#read.csv("16s_soil_otutable_rarefied_nosingleton_fa19_ajo.shared.botu", sep='\t') %>% dim()
# To confirm that numOtus really corresponds to number of Otu column names
#read.csv("16s.bark.sub.shared.botu", sep='\t') %>% dim() # 3282 - 3 columns (label, Group, numBotus) = 3279

# Read in bacterial shaving data and load it as a phyloseq object
bshavings.names <- (read_tsv("16s_bark_otutable_rarefied_nosingleton_fa19_ajo.shared.botu", col_types=cols(Group=col_character()) ))$Group
bshavings.ps <- import_mothur(mothur_shared_file = "16s_bark_otutable_rarefied_nosingleton_fa19_ajo.shared.botu", mothur_constaxonomy_file = "16s_bark_taxonomy_fa19_AJO.cons.taxonomy.botu")

# Reformat sample names so that they are consistent between fungi and bacteria and to match metadata file
sample_names(fshavings.ps) <- fshavings.names %>% (function (x) gsub("-", "_", x)) %>% (function (x) gsub("_L001", "", x)) %>% (function (x) gsub("_001", "", x)) %>% (function (x) gsub("_R1", "", x)) #%>% (function (x) gsub("_S[0-9]{1,2}", "", x))
sample_names(bshavings.ps) <- sapply(bshavings.names, function (x) get_metadata_bacteria()$sample[which(get_metadata_bacteria()$Bact.names == x)])

# Subset the metadata file to match rarified dataset if samples were dropped
fsampdat <- sample_data(inner_join(get_metadata(), data.frame(sample = sample_names(fshavings.ps))))
rownames(fsampdat) <- fsampdat$sample
bsampdat <- sample_data(inner_join(get_metadata_bacteria(), data.frame(sample = sample_names(bshavings.ps))))
rownames(bsampdat) <- bsampdat$sample

# Create a list of consensus samples names that were kept from bacterial and fungal rarefactions
consensus.sampnames.shavings <- data.frame(sample = intersect(fsampdat$sample,bsampdat$sample))

# Shavings.merge is an inner join between shared file, taxonomy file, and metadata file
# cbinded to bacteria shared, taxonomy, and metadata files
fshavings.merge <- merge_phyloseq(fshavings.ps, fsampdat[fsampdat$sample %in% consensus.sampnames.shavings$sample,])
bshavings.merge <- merge_phyloseq(bshavings.ps, bsampdat[bsampdat$sample %in% consensus.sampnames.shavings$sample,])
shavings.merge <- merge_phyloseq(fshavings.merge, bshavings.merge)

# These taxonomy names are not helpful, so let's rename them
colnames(tax_table(shavings.merge)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

## Now repeat for soil

fsoils.names <- (read_tsv("Soil.ITS2.jn.Oct2019.opti_mcc.0.02.subsample.shared", col_types=cols(Group=col_character())) )$Group %>% sub('_S_ITS.+', '', ., perl=T)
fsoils.ps <- import_mothur(mothur_shared_file = "Soil.ITS2.jn.Oct2019.opti_mcc.0.02.subsample.shared", mothur_constaxonomy_file = "soil.mb.its2.good.unique.precluster.pick.pick.opti_mcc.0.02.cons.taxonomy")

bsoils.names <- (read_tsv("16s_soil_otutable_rarefied_nosingleton_fa19_ajo.shared.botu", col_types=cols(Group=col_character()) ))$Group
bsoils.ps <- import_mothur(mothur_shared_file = "16s_soil_otutable_rarefied_nosingleton_fa19_ajo.shared.botu", mothur_constaxonomy_file = "16s_soil_taxonomy_fa19_AJO.cons.taxonomy.botu")
sample_names(fsoils.ps) <- fsoils.names %>% (function (x) gsub("-", "_", x)) %>% (function (x) gsub("_L001", "", x)) %>% (function (x) gsub("_001", "", x)) %>% (function (x) gsub("_R1", "", x))
sample_names(bsoils.ps) <- sapply(bsoils.names, function (x) get_metadata_bacteria()$New.names[which(get_metadata_bacteria()$Bact.names.2 == x)])
get_metadata_soilmb()$sample
get_metadata_bacteria()$New.names
fsampdat <- sample_data(inner_join(get_metadata_soilmb(), data.frame(sample = sample_names(fsoils.ps))))
rownames(fsampdat) <- fsampdat$sample
bsampdat <- sample_data(inner_join(get_metadata_bacteria(), data.frame(New.names = sample_names(bsoils.ps))))
rownames(bsampdat) <- bsampdat$New.names
bsampdat$sample <- bsampdat$New.names
consensus.sampnames.soil <- data.frame(sample = intersect(fsampdat$sample,bsampdat$sample))
fsoils.merge <- merge_phyloseq(fsoils.ps, fsampdat[fsampdat$sample %in% consensus.sampnames.soil$sample,])
bsoils.merge <- merge_phyloseq(bsoils.ps, bsampdat[bsampdat$sample %in% consensus.sampnames.soil$sample,])
soils.merge <- merge_phyloseq(fsoils.merge, bsoils.merge)
colnames(tax_table(soils.merge)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

## Now save global environment to be used in future scripts so we don't have to run this script ever again
## Saves phyloseq object with needed mothur data loaded into the R Environment
setwd("..")
save.image("R_Environments/Jnigra.microbiome.merged.WilliamsOnufrak.RData")
