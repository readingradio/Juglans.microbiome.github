setwd(".")

library(tidyverse)
library(Hmisc)
library(MASS)
library(igraph)
library(tidygraph)
library(ggraph)
library(gridExtra)
library(phyloseq)

source("Net.complex.functions.R")

load("R_Environments/Jnigra.microbiome.merged.WilliamsOnufrak.RData")

soils.merge.df <- otu_table(soils.merge) %>% t() %>% as.data.frame()

state.f <- regmatches(row.names(soils.merge.df), regexpr("(IN|TN|WA)", row.names(soils.merge.df), perl=T))

for (cur.s in c("WA","IN","TN")) {

  ## select only those OTUS that occur in WA/IN/TN at least 5 times in at least 5 samples

  xx <- sapply(1:dim(soils.merge.df)[2], FUN= function (x) sum( soils.merge.df[state.f==cur.s, x] > 5 ) > 5)

  ###
  ### make correlation table with merged soils and soils data
  
  r_otu.all <- rcorr(as.matrix(soils.merge.df[state.f==cur.s,xx]), type='spearman')
  
  triangle <- upper.tri(r_otu.all$r,diag=T)
  r_otu.all$r[triangle] <- NA
  r_otu.all$P[triangle] <- NA

  r.otus <-reshape2::melt(r_otu.all$r, na.rm=T, value.name="r")
  p.otus <-reshape2::melt(r_otu.all$P, na.rm=T, value.name="p")

  cyto.unlisted.all <- cbind(r.otus, data.frame(p=p.otus$p))

  cyto.unlisted.all$r.sign <- as.numeric(cyto.unlisted.all$r > 0)

  setwd(paste("./Net.analysis/CorrNetworks.", cur.s, sep=""))

  # correlation networks
  
  cyto.unlisted.all[cyto.unlisted.all$p < 0.00001 & abs(cyto.unlisted.all$r) > 0.6,]  %>% write.table(., 'cyto.soils.merged.p.0.00001.r.0.6.txt', quote=F, sep="	", row.names=F)
  cyto.unlisted.all[cyto.unlisted.all$p < 0.00001 & abs(cyto.unlisted.all$r) > 0.8,]  %>% write.table(., 'cyto.soils.merged.p.0.00001.r.0.8.txt', quote=F, sep="	", row.names=F)
  cyto.unlisted.all[cyto.unlisted.all$p < 0.00005 & abs(cyto.unlisted.all$r) > 0.6,]  %>% write.table(., 'cyto.soils.merged.p.0.00005.r.0.6.txt', quote=F, sep="	", row.names=F)
  cyto.unlisted.all[cyto.unlisted.all$p < 0.00005 & abs(cyto.unlisted.all$r) > 0.8,]  %>% write.table(., 'cyto.soils.merged.p.0.00005.r.0.8.txt', quote=F, sep="	", row.names=F)
  cyto.unlisted.all[cyto.unlisted.all$p < 0.000005 & abs(cyto.unlisted.all$r) > 0.6,] %>% write.table(., 'cyto.soils.merged.p.0.000005.r.0.6.txt', quote=F, sep="	", row.names=F)
  cyto.unlisted.all[cyto.unlisted.all$p < 0.000005 & abs(cyto.unlisted.all$r) > 0.8,] %>% write.table(., 'cyto.soils.merged.p.0.000005.r.0.8.txt', quote=F, sep="	", row.names=F)
  cyto.unlisted.all[cyto.unlisted.all$p < 0.000001 & abs(cyto.unlisted.all$r) > 0.6,] %>% write.table(., 'cyto.soils.merged.p.0.000001.r.0.6.txt', quote=F, sep="	", row.names=F)
  cyto.unlisted.all[cyto.unlisted.all$p < 0.000001 & abs(cyto.unlisted.all$r) > 0.8,] %>% write.table(., 'cyto.soils.merged.p.0.000001.r.0.8.txt', quote=F, sep="	", row.names=F)

  data.bd <- data.frame(name=NULL, BetweennessCentrality=NULL, Degree=NULL, p=NULL, r=NULL)

  for (i in c( "0.00005","0.00001","0.000005","0.000001")) {
  for (cutoff in c("0.6","0.8")) {
    
	  gr <-read.csv(paste("cyto.soils.merged.p.", i, ".r.", cutoff,".txt", sep = ''), sep='\t') %>% as.matrix() %>% graph_from_data_frame(directed =F)
	  cyto <- data.frame(Degree = degree(gr))
	  cyto$BetweennessCentrality = betweenness(gr)

	  degree.w<-fitdistr(cyto$Degree, 'Weibull')
	  cutoff.d<-qweibull(0.9, degree.w$estimate[1], degree.w$estimate[2])
	  out.d <- cyto$Degree > cutoff.d

  	betweenness.exp <- fitdistr(cyto$BetweennessCentrality, "exponential")
  	cutoff.b<-qexp(0.9, rate = betweenness.exp $estimate)
	  out.b <- cyto$BetweennessCentrality > cutoff.b
	
  	prefix <- paste("cyto.soils.merged.p.", i, ".r.", cutoff, ".network.analysis", sep = '')

	  cyto[out.b & out.d,] %>% write.csv(.,paste(prefix, '.bd.csv', sep=''), quote=F,row.names=T)
	  cyto$name=rownames(cyto)
	
	  data.bd <- rbind(data.bd, cbind(cyto[out.d & out.b,c('name','BetweennessCentrality','Degree')], p=i, r=cutoff), make.row.names=F)
  
	}}

  bd.tab <- with(data.bd, data.frame(name=unique(name)))

  bd.tab$p0.00005.r0.6 <- as.integer(bd.tab$name %in% data.bd[data.bd$p=='0.00005' & data.bd$r=='0.6','name'])
  bd.tab$p0.00001.r0.6 <- as.integer(bd.tab$name %in% data.bd[data.bd$p=='0.00001' & data.bd$r=='0.6','name'])
  bd.tab$p0.000005.r0.6 <- as.integer(bd.tab$name %in% data.bd[data.bd$p=='0.000005' & data.bd$r=='0.6','name'])
  bd.tab$p0.000001.r0.6 <- as.integer(bd.tab$name %in% data.bd[data.bd$p=='0.000001' & data.bd$r=='0.6','name'])

  bd.tab$p0.00005.r0.8 <- as.integer(bd.tab$name %in% data.bd[data.bd$p=='0.00005' & data.bd$r=='0.8','name'])
  bd.tab$p0.00001.r0.8 <- as.integer(bd.tab$name %in% data.bd[data.bd$p=='0.00001' & data.bd$r=='0.8','name'])
  bd.tab$p0.000005.r0.8 <- as.integer(bd.tab$name %in% data.bd[data.bd$p=='0.000005' & data.bd$r=='0.8','name'])
  bd.tab$p0.000001.r0.8 <- as.integer(bd.tab$name %in% data.bd[data.bd$p=='0.000001' & data.bd$r=='0.8','name'])

  bd.tab %>% cbind( tax_table(soils.merge)[as.character(bd.tab$name),], .) %>% write.csv(., paste("../../Results/Soils.HubOTUS.", cur.s, ".csv", sep=""), quote=F, row.names=T)
  
  # matches w top 40 RF OTUs
  bd.tab <- read.csv(paste("../../Results/Soils.HubOTUS.", cur.s, ".csv", sep=""))
  otunames.rf <- read.csv("../../Results/RF.soil.t.table.csv", row.names=1)
  x<- intersect(as.character(bd.tab$name), rownames(otunames.rf))
  otunames.rf[x,] %>% write.csv(., paste("../../Results/Soils.HubOTUS.RF.consensus.", cur.s, ".csv", sep=""), quote=F, row.names=T)

  setwd("../..")

}

## calculate network complexity for each network at each level of p and r

out <-NULL

for (p in c( "0.00005", "0.00001","0.000005","0.000001")) {
  for (r in c("0.6","0.8")) {
  	
  	IN <- read.network (p, r, "IN", subs=T, tissue="soils", path = ".")
	TN <- read.network (p, r, "TN", subs=T, tissue="soils", path = ".")
  	WA <- read.network (p, r, "WA", subs=T, tissue="soils", path = ".")
  	
  	states <- list(IN,TN,WA)
  	names(states) <- c("IN", "TN", "WA")
  	out <- rbind(out, cbind(data.frame(P = rep(p, 3), R = rep(r, 3)), compare.networks(states)))
}}

out$P <- out$P %>% as.character() %>% as.numeric()

out %>% write.csv("Results/Network.complexity.soil.csv")

save.image("R_Environments/Soil.net.stats.RData")

#out <- read.csv("Results/Network.complexity.soil.csv")[,-1]

quartz()
net.complex.multipanel (out, which.leg="first")
net.complex.multipanel (out[out$P<0.00005,], which.leg="first")

## hubs shared between states

in.hubs <- read.csv("Results/Soils.HubOTUS.IN.csv")
tn.hubs <- read.csv("Results/Soils.HubOTUS.TN.csv")
wa.hubs <- read.csv("Results/Soils.HubOTUS.WA.csv")

intersect(in.hubs$name, tn.hubs$name)
intersect(wa.hubs$name, tn.hubs$name)
intersect(wa.hubs$name, in.hubs$name)
intersect(wa.hubs$name, in.hubs$name) %>% intersect(tn.hubs$name)

library(VennDiagram)

venn.diagram(x = list(in.hubs$name, wa.hubs$name, tn.hubs$name), category.names=c("IN","WA","TN"), filename="Figures/Soil.Net.Venn.png", output=T, cat.pos=c(0,0,0), cat.col =tanagr_palette("cyanerpes_cyaneus")[2:4], col=tanagr_palette("cyanerpes_cyaneus")[2:4])
 
## figures
for (p in c( "0.00005", "0.00001","0.000005","0.000001")) {
  for (r in c("0.6","0.8")) {
    
    pdf(paste("Net.analysis/Net.figures/Soil/Allstates.networks.soil.", "p", p, ".r", r, ".nosubnets.pdf", sep=""), width=12, height=6)

    grid.arrange(
      figure.network(p, r, "IN", tissue="soils"),
      figure.network(p, r, "TN", tissue="soils"),
      figure.network(p, r, "WA", tissue="soils"),
      nrow=1,
      top=paste("p = ", p, "r = ", r, sep="")
    )

    dev.off()
    
    pdf(paste("Net.analysis/Net.figures/Soil/Subnetworks/Allstates.networks.soil.", "p", p, ".r", r, ".wsubnets.pdf", sep=""), width=12, height=6)
    
    grid.arrange(
      figure.network(p, r, "IN", T, tissue="soils"),
      figure.network(p, r, "TN", T, tissue="soils"),
      figure.network(p, r, "WA", T, tissue="soils"),
      nrow=1,
      top=paste("p = ", p, "r = ", r, sep="")
    )
    
    dev.off()
  }}

###