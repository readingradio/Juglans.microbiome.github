setwd("/Users/will1809/OneDrive - purdue.edu/TCD microbiome spring 2017/Aaron manuscript/Geoff_Final_RFiles/Juglans.microbiome.github/Mothur_output")

library(tidyverse)
library(readxl)
library(gridExtra)
library(ggpubr)

source('../baxter.nov2019.R')

rarefaction.curve <- function (curves, legpos='none') {
	curves$k <- with(curves, numsampled/1000)
	ggplot(curves, aes(x=k, y=sobs, group=sample, color=State)) +
		geom_line() + 
		theme_bw() +
		theme(
			panel.border = element_blank(),
			panel.background = element_rect(fill = "transparent"), # bg of the panel
    		plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
      	  	legend.position = legpos,
			panel.grid = element_blank(),
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			axis.line = element_line(colour = "black"),
			axis.title = element_text(size = 16),
			axis.text = element_text(size = 14))
}

rare.xlab <- function (main="") labs(title=main, x="\n1000 Sequences Sampled", y="")
rare.ylab <- function (main="") labs(title=main, x="", y="OTU Richness\n")
rare.main <- function (main="") labs(title=main, x="", y="")
rare.corn <- function (main="") labs(title=main, x="\n1000 Sequences Sampled", y="OTU Richness\n")

rare.lgnd <- function () theme(legend.text=element_text(size=20), legend.title=element_text(size=20))

##### ITS2 SOIL

metadata.soil.its2 <- get_metadata_soilmb()

groups.soil.its2  <- read_tsv('Soil.ITS2.jn.Oct2019.opti_mcc.groups.rarefaction')  %>% select(-contains("lci-"), -contains("hci-")) %>% gather(-numsampled, key=sample, value=sobs) %>% mutate(sample=str_replace_all(sample, pattern="0.02-", replacement=""))  %>% drop_na()

# 'sample' has to match

groups.soil.its2$sample[!(groups.soil.its2$sample == 'NTC_ITS_S96')] <- groups.soil.its2$sample[!(groups.soil.its2$sample == 'NTC_ITS_S96')] %>% sub('_S_ITS.+', '', ., perl=T)

# s[!(s == 'NTC_ITS_S96')] %>% sub('_S_ITS.+', '', ., perl=T)

# note drop_na() resulted in 10 fold reduction in number of rows

# change organization of data frames so each observation is a row

its2.soil.curves <- inner_join(metadata.soil.its2, groups.soil.its2)

missing.sample.groups.soil.its2 <- setdiff(groups.soil.its2$sample, metadata.soil.its2$sample)
missing.sample.meta <- setdiff(metadata.soil.its2$sample, groups.soil.its2$sample)
missing.sample.groups.soil.its2
missing.sample.meta

# one of the Washington trees is missing WA132_RN_9


#### ITS2 BRANCHES

metadata.branches.its2 <- get_metadata()

groups.branch.its2 <- read_tsv('Phyllo.ITS2.jn.Oct2019.opti_mcc.branch.groups.rarefaction')  %>% select(-contains("lci-"), -contains("hci-")) %>% gather(-numsampled, key=sample, value=sobs) %>% mutate(sample=str_replace_all(sample, pattern="0.02-", replacement=""))  %>% drop_na()

groups.branch.its2$sample <- gsub("-", "_", groups.branch.its2$sample) %>% (function (x) gsub("_L001", "", x)) %>% (function (x) gsub("_001", "", x)) %>% (function (x) gsub("_R1", "", x))

its2.branch.curves <- inner_join(metadata.branches.its2, groups.branch.its2)

setdiff(groups.branch.its2$sample, its2.branch.curves$sample)
setdiff(alpha.branch.its2$sample, its2.branch.curves$sample)

##### 16S

metadata.bacteria <- get_metadata_bacteria()

##### 16S SOIL

groups.soil.16s  <- read_tsv('soil16sotutable_jnigra17_fa19_ajo.groups.rarefaction')  %>% select(-contains("lci-"), -contains("hci-")) %>% gather(-numsampled, key=sample, value=sobs) %>% mutate(sample=str_replace_all(sample, pattern="0.03-", replacement=""))  %>% drop_na()

metadata.bacteria$sample <- metadata.bacteria$Bact.names.2

B16s.soil.curves <- inner_join(metadata.bacteria, groups.soil.16s )

##### 16S BRANCHES

groups.branches.16s  <- read_tsv('bark16sotutable_jnigra17_fa19_ajo.groups.rarefaction')  %>% select(-contains("lci-"), -contains("hci-")) %>% gather(-numsampled, key=sample, value=sobs) %>% mutate(sample=str_replace_all(sample, pattern="0.03-", replacement=""))  %>% drop_na()

metadata.bacteria$sample <- metadata.bacteria$Bact.names

B16s.branches.curves <- inner_join(metadata.bacteria, groups.branches.16s )

#### PLOTS

its2.soil <- rarefaction.curve(its2.soil.curves)
its2.branches <- rarefaction.curve(its2.branch.curves)
B16s.soil <- rarefaction.curve(B16s.soil.curves)
B16s.branches <- rarefaction.curve(B16s.branches.curves, legpos="none")

grid.arrange(its2.branches + rare.ylab("ITS2 Branches"), its2.soil + rare.main("ITS2 Soil"), B16s.branches + rare.corn("16S Branches"), B16s.soil + rare.xlab("16S Soil"), nrow=2, ncol=2)

plots<-ggarrange(its2.branches + rare.ylab("ITS2 Branches"), its2.soil + rare.main("ITS2 Soil"), B16s.branches + rare.corn("16S Branches"), B16s.soil + rare.xlab("16S Soil"), nrow=2, ncol=2, common.legend= TRUE, legend="right")
ggexport(plots,filename="../Figures/Figure S1.pdf", height=8, width=9)

combined.curves <-
	rbind(data.frame(its2.soil.curves[,c("numsampled","sobs","sample","State")], Compart="Soil", Barcode="ITS2"),
	data.frame(its2.branch.curves[,c("numsampled","sobs","sample","State")], Compart="Branch", Barcode="ITS2"),
	data.frame(B16s.soil.curves[,c("numsampled","sobs","sample","State")], Compart="Soil", Barcode="16S"),
	data.frame(B16s.branches.curves[,c("numsampled","sobs","sample","State")], Compart="Branch", Barcode="16S"))

rarefaction.curve(combined.curves) + facet_grid(Compart ~ Barcode, scales='free')