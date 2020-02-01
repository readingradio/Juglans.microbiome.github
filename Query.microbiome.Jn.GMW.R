# The purpose of this script is to provide tools to explore OTU relationships and abundance of individual taxa in the dataset

library(tidyverse)
library(phyloseq)
library(vegan)

# This function queries a subset of the merged bacterial/fungal dataset by a chosen
# taxonomic level (taxlevel) and the name of the taxon of interest (taxstring).
# It also takes a phyloseq object (pseq).

custom.otu.tax.tbl <- function(pseq, taxlevel, taxstring) {

  onames <- tax_table(pseq)[,taxlevel] %in% taxstring

  rs <- rownames(tax_table(pseq))[onames]

  if (length(rs) < 1) {
    return(NULL)
  } else {
    mtable <- otu_table(pseq)[onames] %>% cbind (tax_table(pseq)[onames], .) %>% t()
    return(mtable)
  }
}

# This function queries the network analysis results for a specified p and r
# cutoff for which a correlation table was generated and gives the taxonomic
# identity of OTUs with significant spearman correlations to G. morbida.

gm.correlations <- function(pseq, pcut, rcut) {
  correlation.network <- read.csv(paste("Net.analysis/CorrNetworks.WA/cyto.shavings.merged.p.", pcut, ".r.", rcut,".txt", sep = ''), sep='\t')
  head(correlation.network)
  v2 <- correlation.network[correlation.network$Var1 == "Otu0115",c('Var2','r','p')]
  v1 <- correlation.network[correlation.network$Var2 == "Otu0115",c('Var1','r','p')]
  return(rbind(tax_table(pseq)[as.character(v2$Var2)] %>% cbind (v2[,2:3]), tax_table(pseq)[as.character(v1$Var1)] %>% cbind (v1[,2:3])))
}

load("R_Environments/Jnigra.microbiome.merged.WilliamsOnufrak.RData")

#Correlated OTUs to G. morbida in network analysis

gm.correlations(shavings.merge, 0.005, 0.6)
gm.correlations(shavings.merge, 0.05, 0.8)
gm.correlations(shavings.merge, 0.01, 0.8)
gm.correlations(shavings.merge, 0.005, 0.8)

# Analysis of relationship between Gm sequence abundance and some other taxa of interest

rownames(custom.otu.tax.tbl(shavings.merge, 'Genus', 'g__Geosmithia'))

gm <- custom.otu.tax.tbl(shavings.merge, 'Genus', 'g__Geosmithia')[-c(1:7),1]
ap <- custom.otu.tax.tbl(shavings.merge, 'Genus', 'g__Aureobasidium')[-c(1:7),1]
sp <- custom.otu.tax.tbl(shavings.merge, 'Genus', 'g__Sydowia')[-c(1:7),1]
xm <- custom.otu.tax.tbl(shavings.merge, 'Family', 'Xanthomonadaceae')[-c(1:7),1]
pm <- custom.otu.tax.tbl(shavings.merge, 'Genus', 'Pseudomonas')[-c(1:7),2] ## Botu0021

shan <- shavings.merge %>% otu_table() %>% t() %>% diversity(., 'shannon')
rich <- shavings.merge %>% otu_table() %>% t() %>% specnumber()
wa <- data.frame(syd = as.numeric(sp[30:39]), a = as.numeric(ap[30:39]), g = as.numeric(gm[30:39]), r = as.numeric(rich[30:39]), s = as.numeric(shan[30:39]), x =as.numeric(xm[30:39]), p =as.numeric(pm[30:39]))


# Linear models
wa.lm.r <- with(wa, lm(r ~ log1p(g))) %>% summary() ## richness was significant
wa.lm.s <- with(wa, lm(s ~ log1p(g))) %>% summary() ## Shannon diversity was significant
wa.lm.a <- with(wa, lm(log1p(a) ~ log1p(g))) %>% summary() ## Auerobasidium BOTU0057 was significant
wa.lm.syd <- with(wa, lm(log1p(syd) ~ log1p(g))) %>% summary()
wa.lm.x <- with(wa, lm(log1p(x) ~ log1p(g))) %>% summary() ## Xanthomonas BOTU0007 was significant
wa.lm.p <- with(wa, lm(log1p(p) ~ log1p(g))) %>% summary()

# Check assumptions
with(wa, lm(r ~ log1p(g))) %>% plot()
with(wa, lm(s ~ log1p(g))) %>% plot()
with(wa, lm(log1p(a) ~ log1p(g))) %>% plot()
with(wa, lm(log1p(x) ~ log1p(g))) %>% plot()

# Final plots

quartz()
par(mfrow=c(1,2))
with(wa, plot(log1p(g),r, xlab = expression(paste("Log reads ", italic("Geosmithia "), "OTU0115 (WA only)")), ylab= "Richness", pch=3))
points(wa$g[order(wa$g)], (wa.lm.r$coefficients[1] + wa.lm.r$coefficients[2] * wa$g)[order(wa$g)], type='l')
text(0,285,paste("F =",wa.lm.r$fstatistic[[1]] %>% round(2)), pos=4)
text(0,280,paste("p =",wa.lm.r$coefficients[[2,4]] %>% round(3)), pos=4)
text(0,275,paste("R² =",wa.lm.r$r.squared %>% round(3)), pos=4)
with(wa, plot(log1p(g),s, xlab = expression(paste("Log reads ", italic("Geosmithia "), "OTU0115 (WA only)")), ylab= "Shannon Diversity", pch=3))
points(wa$g[order(wa$g)], (wa.lm.s$coefficients[1] + wa.lm.s$coefficients[2] * wa$g)[order(wa$g)], type='l')
text(0,4.05,paste("F =",wa.lm.s$fstatistic[[1]] %>% round(2)), pos=4)
text(0,4,paste("p =",wa.lm.s$coefficients[[2,4]] %>% round(3)), pos=4)
text(0,3.95,paste("R² =",wa.lm.s$r.squared %>% round(3)), pos=4)

quartz()
par(mfrow=c(2,2))
with(wa, plot(log1p(g),log1p(a), xlab = expression(paste("Log reads ", italic("Geosmithia "), "OTU0115 (WA only)")), ylab= expression(paste("Log reads ", italic("A. pullulans "), "OTU0057")), pch=3))
points(wa$g[order(wa$g)], (wa.lm.a$coefficients[1] + wa.lm.a$coefficients[2] * wa$g)[order(wa$g)], type='l')
text(2,5.8,paste("F =",wa.lm.a$fstatistic[[1]] %>% round(2)), pos=4)
text(2,5.6,paste("p =",wa.lm.a$coefficients[[2,4]] %>% round(3)), pos=4)
text(2,5.4,paste("R² =",wa.lm.a$r.squared %>% round(3)), pos=4)
with(wa, plot(log1p(g),log1p(syd), xlab = expression(paste("Log reads ", italic("Geosmithia "), "OTU0115 (WA only)")), ylab= expression(paste("Log reads ", italic("Sydowia "), "OTU0194")), pch=3))
text(1.5,2.25,paste("F =",wa.lm.syd$fstatistic[[1]] %>% round(2)), pos=4)
text(1.5,2,paste("p =",wa.lm.syd$coefficients[[2,4]] %>% round(3)), pos=4)
text(1.5,1.75,paste("R² =",wa.lm.syd$r.squared %>% round(3)), pos=4)
with(wa, plot(log1p(g),log1p(x), xlab = expression(paste("Log reads ", italic("Geosmithia "), "OTU0115 (WA only)")), ylab=expression(paste("Log reads ", italic("Xanthomonas "), "BOTU0007")), pch=3))
points(wa$g[order(wa$g)], (wa.lm.x$coefficients[1] + wa.lm.x$coefficients[2] * wa$g)[order(wa$g)], type='l')
text(2,2,paste("F =",wa.lm.x$fstatistic[[1]] %>% round(2)), pos=4)
text(2,1.5,paste("p =",wa.lm.x$coefficients[[2,4]] %>% round(3)), pos=4)
text(2,1,paste("R² =",wa.lm.x$r.squared %>% round(3)), pos=4)
with(wa, plot(log1p(g),log1p(p), xlab = expression(paste("Log reads ", italic("Geosmithia "), "OTU0115 (WA only)")), ylab= expression(paste("Log reads ", italic("Pseudomonas "), "BOTU0021")), pch=3))
text(2,2,paste("F =",wa.lm.p$fstatistic[[1]] %>% round(2)), pos=4)
text(2,1.5,paste("p =",wa.lm.p$coefficients[[2,4]] %>% round(3)), pos=4)
text(2,1,paste("R² =",wa.lm.p$r.squared %>% round(3)), pos=4)
# Look at mean + SD of Gm reads

custom.otu.tax.tbl(shavings.merge, 'Genus', 'g__Geosmithia')
rownames(custom.otu.tax.tbl(shavings.merge, 'Genus', 'g__Geosmithia'))
custom.otu.tax.tbl(shavings.merge, 'Genus', 'g__Geosmithia')[37:46,1] %>% as.numeric() %>% sd()
custom.otu.tax.tbl(shavings.merge, 'Genus', 'g__Geosmithia')[37:46,1] %>% as.numeric() %>% mean()

# Queries

custom.otu.tax.tbl(shavings.merge, 'Genus', 'g__Armillaria')
custom.otu.tax.tbl(shavings.merge, 'Genus', 'g__Phellinus')
custom.otu.tax.tbl(shavings.merge, 'Family', 'f__Nectriaceae')[,1:10]
custom.otu.tax.tbl(shavings.merge, 'Family', 'f__Nectriaceae')[,11:21]
custom.otu.tax.tbl(shavings.merge, 'Genus', 'Pseudomonas')
custom.otu.tax.tbl(shavings.merge, 'Genus', 'Brenneria')
custom.otu.tax.tbl(shavings.merge, 'Genus', 'Erwinia')
custom.otu.tax.tbl(shavings.merge, 'Genus', 'Pectobacterium')
custom.otu.tax.tbl(shavings.merge, 'Genus', 'Agrobacterium')
custom.otu.tax.tbl(shavings.merge, 'Family','Xanthomonadaceae')
custom.otu.tax.tbl(shavings.merge, 'Genus', 'g__Schizophyllum')
custom.otu.tax.tbl(shavings.merge, 'Genus', 'g__Hypochnicium')
custom.otu.tax.tbl(shavings.merge, 'Genus', 'g__Trametes')
custom.otu.tax.tbl(shavings.merge, 'Genus', 'g__Peniophora')
custom.otu.tax.tbl(shavings.merge, 'Genus', 'g__Hericium')

custom.otu.tax.tbl(shavings.merge, 'Genus', 'g__Sydowia')				# One OTU0194 only present in WA (polyspora)
custom.otu.tax.tbl(shavings.merge, 'Genus', 'g__Endoconidioma')			# much more abundant in WA (E. populi OTU0010) but also present in TN * also in RF
custom.otu.tax.tbl(shavings.merge, 'Order', 'o__Botryosphaeriales')		# more abundant in IN and TN (OTU0023 and OTU0065)

custom.otu.tax.tbl(shavings.merge, 'Order', 'o__Capnodiales') %>% dim()
custom.otu.tax.tbl(shavings.merge, 'Order', 'o__Capnodiales')[,1:10]
custom.otu.tax.tbl(shavings.merge, 'Order', 'o__Capnodiales')[,11:20]			# mostly present in TN, lesser extent IN, but OTU0012 Mycosphaerella tassiana much more abundant in WA
custom.otu.tax.tbl(shavings.merge, 'Order', 'o__Capnodiales')[,21:30]
custom.otu.tax.tbl(shavings.merge, 'Order', 'o__Capnodiales')[,31:40]

custom.otu.tax.tbl(shavings.merge, 'Genus', 'g__Aureobasidium') %>% dim()
custom.otu.tax.tbl(shavings.merge, 'Genus', 'g__Aureobasidium') # two OTUs of pullulans much more abundant in WA but one also present at low levels in INTN OTU0006 and one not OTU 0057

custom.otu.tax.tbl(shavings.merge, 'Family', 'f__Didymellaceae')			# OTU0001 much more abundant in WA but also present in IN and TN (genus unk)
custom.otu.tax.tbl(shavings.merge, 'Family', 'f__Lophiostomataceae')	# More abundant in TN (OTU 0101 fuckelii) and IN (OTU0125)
custom.otu.tax.tbl(shavings.merge, 'Family', 'f__Didymellaceae')
custom.otu.tax.tbl(shavings.merge, 'Family', 'f__Lophiostomataceae')
custom.otu.tax.tbl(shavings.merge, 'Family', 'f__Phaeosphaeriaceae')
custom.otu.tax.tbl(shavings.merge, 'Phylum', 'p__Chytridiomycota')	
custom.otu.tax.tbl(shavings.merge, 'Family', 'f__Herpotrichiellaceae')	
custom.otu.tax.tbl(shavings.merge, 'Genus', 'Pseudomonas')	
custom.otu.tax.tbl(shavings.merge, 'Family', 'Burkholderiaceae')	
custom.otu.tax.tbl(shavings.merge, 'Order', 'Betaproteobacteriales ')	
custom.otu.tax.tbl(shavings.merge, 'Order', 'o__Eurotiales')
custom.otu.tax.tbl(shavings.merge, 'Order', 'o__Chaetothyriales')[,1:5]
custom.otu.tax.tbl(shavings.merge, 'Order', 'o__Chaetothyriales')[,6:10]
custom.otu.tax.tbl(shavings.merge, 'Order', 'o__Phaeomoniellales')[,1:5]
custom.otu.tax.tbl(shavings.merge, 'Order', 'o__Phaeomoniellales')[,6:10]
custom.otu.tax.tbl(shavings.merge, 'Order', 'o__Xylariales')[,1:5]
custom.otu.tax.tbl(shavings.merge, 'Order', 'o__Xylariales')[,6:10]
custom.otu.tax.tbl(shavings.merge, 'Genus', 'g__Knufia')
custom.otu.tax.tbl(shavings.merge, 'Genus', 'g__Colletotrichum')
custom.otu.tax.tbl(shavings.merge, 'Order', 'o__Hypocreales')[,1:5]
custom.otu.tax.tbl(shavings.merge, 'Order', 'o__Xylariales')[,6:10]
custom.otu.tax.tbl(shavings.merge, 'Class', 'c__Sordariomycetes')[,1:10]
custom.otu.tax.tbl(shavings.merge, 'Class', 'c__Tremellomycetes')[,1:10]
custom.otu.tax.tbl(shavings.merge, 'Class', 'c__Taphrinomycetes')[,11:20]
custom.otu.tax.tbl(shavings.merge, 'Family', 'Chitinophagaceae')
custom.otu.tax.tbl(shavings.merge, 'Order', 'o__Botryosphaeriales')		# more abundant in IN and TN (OTU0023 and OTU0065)
custom.otu.tax.tbl(shavings.merge, 'Order', 'o__Ophiostomatales')

#########################
## SOIL #################
#########################

custom.otu.tax.tbl(soils.merge, 'Genus', 'g__Fusarium')[,1:10]
custom.otu.tax.tbl(soils.merge, 'Family', 'f__Xylariaceae')
custom.otu.tax.tbl(soils.merge, 'Genus', 'g__Neocosmospora')
custom.otu.tax.tbl(soils.merge, 'Genus', 'g__Cylindrocladium')

custom.otu.tax.tbl(soils.merge, 'Genus', 'g__Armillaria')
custom.otu.tax.tbl(soils.merge, 'Genus', 'g__Phellinus')
custom.otu.tax.tbl(soils.merge, 'Genus', 'Pseudomonas')
custom.otu.tax.tbl(soils.merge, 'Genus', 'g__Schizophyllum')
custom.otu.tax.tbl(soils.merge, 'Genus', 'g__Hypochnicium')
custom.otu.tax.tbl(soils.merge, 'Genus', 'g__Trametes')
custom.otu.tax.tbl(soils.merge, 'Genus', 'g__Peniophora')
custom.otu.tax.tbl(soils.merge, 'Genus', 'g__Hericium')
custom.otu.tax.tbl(soils.merge, 'Genus', 'g__Pleurotus')

custom.otu.tax.tbl(soils.merge, 'Family', 'f__Nectriaceae')[,1:20]
custom.otu.tax.tbl(soils.merge, 'Phylum', 'p__Glomeromycota')[,1:20]

custom.otu.tax.tbl(soils.merge, 'Family', 'f__Nectriaceae')[8:52, ]  %>% apply(., 2, as.numeric) %>% diversity() %>% data.frame(State = c(rep("IN", 15), rep("TN", 18), rep("WA", 12)), Diversity=.) %>% lm(Diversity ~ State, .) %>% anova()
custom.otu.tax.tbl(soils.merge, 'Family', 'f__Nectriaceae')[8:52, ]  %>% apply(., 2, as.numeric) %>% diversity() %>% data.frame(State = c(rep("IN", 15), rep("TN", 18), rep("WA", 12)), Diversity=.) %>% aov(Diversity ~ State, .) %>% TukeyHSD()
custom.otu.tax.tbl(soils.merge, 'Family', 'f__Nectriaceae')[8:52, ]  %>% apply(., 2, as.numeric) %>% specnumber() %>% data.frame(State = c(rep("IN", 15), rep("TN", 18), rep("WA", 12)), Diversity=.) %>% lm(Diversity ~ State, .) %>% anova()
custom.otu.tax.tbl(soils.merge, 'Family', 'f__Nectriaceae')[8:52, ]  %>% apply(., 2, as.numeric) %>% specnumber() %>% data.frame(State = c(rep("IN", 15), rep("TN", 18), rep("WA", 12)), Diversity=.) %>% aov(Diversity ~ State, .) %>% TukeyHSD()

custom.otu.tax.tbl(soils.merge, 'Genus', c('g__Fusarium', 'g__Cylindrocarpon'))[8:52, ]  %>% apply(., 2, as.numeric) %>% diversity() %>% data.frame(State = c(rep("IN", 15), rep("TN", 18), rep("WA", 12)), Diversity=.) %>% lm(Diversity ~ State, .) %>% anova()
custom.otu.tax.tbl(soils.merge, 'Genus', c('g__Fusarium', 'g__Cylindrocarpon'))[8:52, ]  %>% apply(., 2, as.numeric) %>% diversity() %>% data.frame(State = c(rep("IN", 15), rep("TN", 18), rep("WA", 12)), Diversity=.) %>% aov(Diversity ~ State, .) %>% TukeyHSD()
custom.otu.tax.tbl(soils.merge, 'Genus', c('g__Fusarium', 'g__Cylindrocarpon'))[8:52, ]  %>% apply(., 2, as.numeric) %>% specnumber() %>% data.frame(State = c(rep("IN", 15), rep("TN", 18), rep("WA", 12)), Diversity=.) %>% lm(Diversity ~ State, .) %>% anova()
custom.otu.tax.tbl(soils.merge, 'Genus', c('g__Fusarium', 'g__Cylindrocarpon'))[8:52, ]  %>% apply(., 2, as.numeric) %>% specnumber() %>% data.frame(State = c(rep("IN", 15), rep("TN", 18), rep("WA", 12)), Diversity=.) %>% aov(Diversity ~ State, .) %>% TukeyHSD()

custom.otu.tax.tbl(soils.merge, 'Phylum', 'p__Glomeromycota')[8:52, ]  %>% apply(., 2, as.numeric) %>% diversity() %>% data.frame(State = c(rep("IN", 15), rep("TN", 18), rep("WA", 12)), Diversity=.) %>% lm(Diversity ~ State, .) %>% anova()
custom.otu.tax.tbl(soils.merge, 'Phylum', 'p__Glomeromycota')[8:52, ]  %>% apply(., 2, as.numeric) %>% diversity() %>% data.frame(State = c(rep("IN", 15), rep("TN", 18), rep("WA", 12)), Diversity=.) %>% aov(Diversity ~ State, .) %>% TukeyHSD()

custom.otu.tax.tbl(soils.merge, 'Phylum', 'p__Glomeromycota')[8:52, ]  %>% apply(., 2, as.numeric) %>% specnumber() %>% data.frame(State = c(rep("IN", 15), rep("TN", 18), rep("WA", 12)), Diversity=.) %>% lm(Diversity ~ State, .) %>% anova()
custom.otu.tax.tbl(soils.merge, 'Phylum', 'p__Glomeromycota')[8:52, ]  %>% apply(., 2, as.numeric) %>% specnumber() %>% data.frame(State = c(rep("IN", 15), rep("TN", 18), rep("WA", 12)), Diversity=.) %>% aov(Diversity ~ State, .) %>% TukeyHSD()

custom.otu.tax.tbl(soils.merge, 'Family', 'f__Nectriaceae')[8:52, ]  %>% apply(., 2, as.numeric) 
custom.otu.tax.tbl(soils.merge, 'Phylum', 'p__Glomeromycota')[,1:20]

nectriaceae.pcoa <- custom.otu.tax.tbl(soils.merge, 'Family', 'f__Nectriaceae')[8:52, ]  %>% apply(., 2, as.numeric) %>% vegdist() %>% pcoa()
custom.otu.tax.tbl(soils.merge, 'Family', 'f__Nectriaceae')[8:52, ]  %>% apply(., 2, as.numeric) %>% vegdist() %>% pcoa()

fus.cylincpn.pcoa <- custom.otu.tax.tbl(soils.merge, 'Genus', c('g__Fusarium', 'g__Cylindrocarpon'))[8:52, ]  %>% apply(., 2, as.numeric) %>% vegdist() %>% pcoa()

plot(wa$g[c(1:5,7:8)], nectriaceae.pcoa$vectors[c(38,43:45,41,35,36),'Axis.1'])

rownames(nectriaceae.pcoa$vectors)<-colnames(otu_table(soils.merge))
rownames(nectriaceae.pcoa$vectors) %>% length()
grps <- data.frame(State= as.factor(c(rep("IN", 15), rep("TN", 18), rep("WA", 12))), row.names = colnames(otu_table(soils.merge)))

quartz()
par(mfrow=c(1,2), oma=c(1,1,1,1), xpd=TRUE, mar=c(4,4,2,2))

nectriaceae.pcoa$vectors[,c("Axis.1","Axis.2")] %>% plot(col=c(rep("green", 15), rep("blue", 18), rep("orange", 12)),pch=15, axes=F, xlab="PCOA 1", ylab="PCOA 2", xlim=c(-.6,.6), ylim=c(-.6,.6))
axis(2, line=0);axis(1, line=0)
nectriaceae.pcoa$scores <- nectriaceae.pcoa$vectors
nectriaceae.pcoa %>% ordiellipse(., groups=grps$State, kind = "sd", draw = "lines", col=c("green","blue", "orange"))
text(-.5,.5,"A")

glomeromycota.pcoa <- custom.otu.tax.tbl(soils.merge, 'Phylum', 'p__Glomeromycota')[8:52, ]  %>% apply(., 2, as.numeric) %>% vegdist() %>% pcoa()
glomeromycota.pcoa$vectors[,c("Axis.1","Axis.2")] %>% plot(col=c(rep("green", 15), rep("blue", 18), rep("orange", 12)),pch=15, axes=F, xlab="PCOA 1", ylab="", xlim=c(-.6,.6), ylim=c(-.6,.6))
axis(1, line=0);axis(2, line=0)
glomeromycota.pcoa$scores<-glomeromycota.pcoa$vectors
glomeromycota.pcoa %>% ordiellipse(., groups=grps$State, kind = "sd", draw = "lines", col=c("green","blue", "orange"))
text(-.5,.5,"B")
legend("topright", legend=c("IN","TN","WA"), pch=15, col=c("green","blue", "orange"),horiz=TRUE, bty="n", xpd=NA, inset=c(0,-.1))

quartz()
par(mfrow=c(1,2), oma=c(1,1,1,1), xpd=TRUE, mar=c(4,4,2,2))

fus.cylincpn.pcoa$vectors[,c("Axis.1","Axis.2")] %>% plot(col=c(rep("green", 15), rep("blue", 18), rep("orange", 12)),pch=15, axes=F, xlab="PCOA 1", ylab="PCOA 2", xlim=c(-.6,.6), ylim=c(-.6,.6))
axis(2, line=0);axis(1, line=0)
fus.cylincpn.pcoa$scores <- fus.cylincpn.pcoa $vectors
fus.cylincpn.pcoa  %>% ordiellipse(., groups=grps$State, kind = "sd", draw = "lines", col=c("green","blue", "orange"))
text(-.5,.5,"A")

glomeromycota.pcoa <- custom.otu.tax.tbl(soils.merge, 'Phylum', 'p__Glomeromycota')[8:52, ]  %>% apply(., 2, as.numeric) %>% vegdist() %>% pcoa()
glomeromycota.pcoa$vectors[,c("Axis.1","Axis.2")] %>% plot(col=c(rep("green", 15), rep("blue", 18), rep("orange", 12)),pch=15, axes=F, xlab="PCOA 1", ylab="", xlim=c(-.6,.6), ylim=c(-.6,.6))
axis(1, line=0);axis(2, line=0)
glomeromycota.pcoa$scores<-glomeromycota.pcoa$vectors
glomeromycota.pcoa %>% ordiellipse(., groups=grps$State, kind = "sd", draw = "lines", col=c("green","blue", "orange"))
text(-.5,.5,"B")
legend("topright", legend=c("IN","TN","WA"), pch=15, col=c("green","blue", "orange"),horiz=TRUE, bty="n", xpd=NA, inset=c(0,-.1))

#plot (wa$g[c(1:5,7:8)], glomeromycota.pcoa$vectors[c(38,43:45,41,35,36),'Axis.1'])

lm (glomeromycota.pcoa$vectors[c(38,43:45,41,35,36),'Axis.4'] ~ wa$g[c(1:5,7:8)]) %>% summary()

State = data.frame(s = c(rep("IN", 15), rep("TN", 18), rep("WA", 12)))
nectriaceae.dist <- custom.otu.tax.tbl(soils.merge, 'Family', 'f__Nectriaceae')[8:52, ]  %>% apply(., 2, as.numeric) %>% vegdist()
glomeromycota.dist <- custom.otu.tax.tbl(soils.merge, 'Phylum', 'p__Glomeromycota')[8:52, ]  %>% apply(., 2, as.numeric) %>% vegdist()
adonis(nectriaceae.dist ~ s, data=State, permutations=9999)
adonis(glomeromycota.dist ~ s, data=State, permutations=9999)
fus.cylincpn.dist <- custom.otu.tax.tbl(soils.merge, 'Genus', c('g__Fusarium', 'g__Cylindrocarpon'))[8:52, ]  %>% apply(., 2, as.numeric) %>% vegdist()
adonis(fus.cylincpn.dist ~ s, data=State, permutations=9999)

