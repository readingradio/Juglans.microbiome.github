setwd("/Users/will1809/OneDrive - purdue.edu/TCD microbiome spring 2017/Aaron manuscript/Geoff_Final_RFiles/Juglans.microbiome.github")

library(tidyverse)
library(Hmisc)
library(MASS)
library(igraph)
library(tidygraph)
library(ggraph)
library(gridExtra)
library(phyloseq)

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
  
  # matches

  otunames.rf <- read.csv("../../Results/RF.soil.t.table.csv", row.names=1)
  otunames.rf[as.character(bd.tab$name),] %>% write.csv(., paste("../../Results/Soils.HubOTUS.RF.consensus.", cur.s, ".csv", sep=""), quote=F, row.names=T)

  setwd("../..")

}



## figures

figure.network <- function (p, r, s, subs=F) {

  gr <- read.csv(paste("Net.analysis/CorrNetworks.", s,"/cyto.soils.merged.p.", p, ".r.", r,".txt", sep = ''), sep='\t')[,1:3] %>% as.matrix() %>% graph_from_data_frame(directed =F)
  edge_attr(gr, name= 'r.sign') <- read.csv(paste("Net.analysis/CorrNetworks.", s,"/cyto.soils.merged.p.", p, ".r.", r,".txt", sep = ''), sep='\t')[,3] %>% sign() %>% as.factor()
  edge_attr(gr, name= 'r.abs') <- read.csv(paste("Net.analysis/CorrNetworks.", s,"/cyto.soils.merged.p.", p, ".r.", r,".txt", sep = ''), sep='\t')[,3] %>% abs() %>% as.numeric()
  vertex_attr(gr, name='deg') <- degree(gr)
  vertex_attr(gr, name='fu.ba') <- (vertex_attr(gr, 'name') %>% grepl('Botu',.) ) %>% as.numeric()
  
  bd.tab <- read.csv(paste("Results/Soils.HubOTUS.", s, ".csv", sep=""), row.names = 1)
  vertex_attr(gr, name='hub') <- vertex_attr(gr, 'name') %in% rownames(bd.tab)[bd.tab[, paste("p",p,".r",r,sep="")] == 1] %>% as.numeric()
  
  if (s == "WA") {
    vertex_attr(gr, name='hub') <- ((((vertex_attr(gr, 'name') == 'Otu0115' )) %>% as.numeric) + vertex_attr(gr, name='hub'))
    vertex_attr(gr, name='fu.ba') <- ((((vertex_attr(gr, 'name') == 'Otu0115' )) %>% as.numeric)*2 + vertex_attr(gr, name='fu.ba'))
  }
  
  if (!subs) {
    dec.subrs <- decompose.graph(gr)
    cl.gr.i <- which.max(clusters(gr)$csize)
    gr <- dec.subrs[[cl.gr.i]]
  }
  
  layout <- layout_with_kk(gr, weights=2 - as.numeric(edge_attr(gr, 'r')))
  
  ggraph(gr, layout) + 
    geom_edge_link(aes(width=.1, alpha=r.abs, color=as.factor(r.sign))) + 
    scale_edge_color_manual(values=c('grey','blue')) + 
    scale_edge_width(range=c(0.2,1)) + 
    geom_node_point(aes(size=deg, shape=as.factor(hub), color=as.factor(fu.ba))) + 
    scale_size(range = c(.2,2.5)) + scale_color_manual(values=c('red','black','green')) +
    scale_shape_manual(values=c(1,19,8)) + 
    theme_graph() + 
    theme(legend.position='none')
  
}

for (p in c( "0.00005", "0.00001","0.000005","0.000001")) {
  for (r in c("0.6","0.8")) {
    
    pdf(paste("Net.analysis/Net.figures/Soil/Allstates.networks.soil.", "p", p, ".r", r, ".nosubnets.pdf", sep=""), width=12, height=6)

    grid.arrange(
      figure.network(p, r, "IN"),
      figure.network(p, r, "TN"),
      figure.network(p, r, "WA"),
      nrow=1,
      top=paste("p = ", p, "r = ", r, sep="")
    )

    dev.off()
    
    pdf(paste("Net.analysis/Net.figures/Soil/Subnetworks/Allstates.networks.soil.", "p", p, ".r", r, ".wsubnets.pdf", sep=""), width=12, height=6)
    
    grid.arrange(
      figure.network(p, r, "IN", T),
      figure.network(p, r, "TN", T),
      figure.network(p, r, "WA", T),
      nrow=1,
      top=paste("p = ", p, "r = ", r, sep="")
    )
    
    dev.off()
  }}

###