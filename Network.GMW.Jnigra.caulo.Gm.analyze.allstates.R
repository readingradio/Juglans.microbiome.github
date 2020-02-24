setwd(".")

library(tidyr)
library(Hmisc)
library(MASS)
library(igraph)
library(tidygraph)
library(ggraph)
library(gridExtra)
library(QuACN)
library(acss)

source("Net.complex.functions.R")

load("R_Environments/Jnigra.microbiome.merged.WilliamsOnufrak.RData")

shavings.merge.df <- otu_table(shavings.merge) %>% t() %>% as.data.frame()

state.f <- regmatches(row.names(shavings.merge.df), regexpr("(IN|TN|WA)", row.names(shavings.merge.df), perl=T))

for (cur.s in c("WA","IN","TN")) {

  ## select only those OTUS that occur in WA/IN/TN at least 1 time in at least 3 samples

  xx <- sapply(1:dim(shavings.merge.df)[2], FUN= function (x) sum( shavings.merge.df[state.f==cur.s, x] > 0 ) > 3)

  ###
  ### make correlation table with merged shavings and shavings data
  
  r_otu.all <- rcorr(as.matrix(shavings.merge.df[state.f==cur.s,xx]), type='spearman')
  
  triangle <- upper.tri(r_otu.all$r,diag=T)
  r_otu.all$r[triangle] <- NA
  r_otu.all$P[triangle] <- NA

  r.otus <-reshape2::melt(r_otu.all$r, na.rm=T, value.name="r")
  p.otus <-reshape2::melt(r_otu.all$P, na.rm=T, value.name="p")

  cyto.unlisted.all <- cbind(r.otus, data.frame(p=p.otus$p))

  cyto.unlisted.all$r.sign <- as.numeric(cyto.unlisted.all$r > 0)

  setwd(paste("./Net.analysis/CorrNetworks.", cur.s, sep=""))

  # correlation networks
  cyto.unlisted.all[cyto.unlisted.all$p < 0.05 & abs(cyto.unlisted.all$r) > 0.6,] %>% write.table(., 'cyto.shavings.merged.p.0.05.r.0.6.txt', quote=F, sep="	", row.names=F)
  cyto.unlisted.all[cyto.unlisted.all$p < 0.05 & abs(cyto.unlisted.all$r) > 0.8,] %>% write.table(., 'cyto.shavings.merged.p.0.05.r.0.8.txt', quote=F, sep="	", row.names=F)
  cyto.unlisted.all[cyto.unlisted.all$p < 0.01 & abs(cyto.unlisted.all$r) > 0.6,] %>% write.table(., 'cyto.shavings.merged.p.0.01.r.0.6.txt', quote=F, sep="	", row.names=F)
  cyto.unlisted.all[cyto.unlisted.all$p < 0.01 & abs(cyto.unlisted.all$r) > 0.8,] %>% write.table(., 'cyto.shavings.merged.p.0.01.r.0.8.txt', quote=F, sep="	", row.names=F)
  cyto.unlisted.all[cyto.unlisted.all$p < 0.005 & abs(cyto.unlisted.all$r) > 0.6,] %>% write.table(., 'cyto.shavings.merged.p.0.005.r.0.6.txt', quote=F, sep="	", row.names=F)
  cyto.unlisted.all[cyto.unlisted.all$p < 0.005 & abs(cyto.unlisted.all$r) > 0.8,] %>% write.table(., 'cyto.shavings.merged.p.0.005.r.0.8.txt', quote=F, sep="	", row.names=F)
  cyto.unlisted.all[cyto.unlisted.all$p < 0.001 & abs(cyto.unlisted.all$r) > 0.6,] %>% write.table(., 'cyto.shavings.merged.p.0.001.r.0.6.txt', quote=F, sep="	", row.names=F)
  cyto.unlisted.all[cyto.unlisted.all$p < 0.001 & abs(cyto.unlisted.all$r) > 0.8,] %>% write.table(., 'cyto.shavings.merged.p.0.001.r.0.8.txt', quote=F, sep="	", row.names=F)

  data.bd <- data.frame(name=NULL, BetweennessCentrality=NULL, Degree=NULL, p=NULL, r=NULL)

  for (i in c("0.05","0.01","0.005", "0.001")) {
  for (cutoff in c("0.6","0.8")) {
    
	  gr <-read.csv(paste("cyto.shavings.merged.p.", i, ".r.", cutoff,".txt", sep = ''), sep='\t') %>% as.matrix() %>% graph_from_data_frame(directed =F)
	  cyto <- data.frame(Degree = degree(gr))
	  cyto$BetweennessCentrality = betweenness(gr)

	  degree.w<-fitdistr(cyto$Degree, 'Weibull')
	  cutoff.d<-qweibull(0.9, degree.w$estimate[1], degree.w$estimate[2])
	  out.d <- cyto$Degree > cutoff.d

  	betweenness.exp <- fitdistr(cyto$BetweennessCentrality, "exponential")
  	cutoff.b<-qexp(0.9, rate = betweenness.exp $estimate)
	  out.b <- cyto$BetweennessCentrality > cutoff.b
	
  	prefix <- paste("cyto.shavings.merged.p.", i, ".r.", cutoff, ".network.analysis", sep = '')

	  cyto[out.b & out.d,] %>% write.csv(.,paste(prefix, '.bd.csv', sep=''), quote=F,row.names=T)
	  cyto$name=rownames(cyto)
	
	  data.bd <- rbind(data.bd, cbind(cyto[out.d & out.b,c('name','BetweennessCentrality','Degree')], p=i, r=cutoff), make.row.names=F)
  
	}}

  bd.tab <- with(data.bd, data.frame(name=unique(name)))

  bd.tab$p0.05.r0.6 <- as.integer(bd.tab$name %in% data.bd[data.bd$p=='0.05' & data.bd$r=='0.6','name'])
  bd.tab$p0.01.r0.6 <- as.integer(bd.tab$name %in% data.bd[data.bd$p=='0.01' & data.bd$r=='0.6','name'])
  bd.tab$p0.005.r0.6 <- as.integer(bd.tab$name %in% data.bd[data.bd$p=='0.005' & data.bd$r=='0.6','name'])
  bd.tab$p0.001.r0.6 <- as.integer(bd.tab$name %in% data.bd[data.bd$p=='0.0001' & data.bd$r=='0.6','name'])

  bd.tab$p0.05.r0.8 <- as.integer(bd.tab$name %in% data.bd[data.bd$p=='0.05' & data.bd$r=='0.8','name'])
  bd.tab$p0.01.r0.8 <- as.integer(bd.tab$name %in% data.bd[data.bd$p=='0.01' & data.bd$r=='0.8','name'])
  bd.tab$p0.005.r0.8 <- as.integer(bd.tab$name %in% data.bd[data.bd$p=='0.005' & data.bd$r=='0.8','name'])
  bd.tab$p0.001.r0.8 <- as.integer(bd.tab$name %in% data.bd[data.bd$p=='0.0001' & data.bd$r=='0.8','name'])

  bd.tab %>% cbind( tax_table(shavings.merge)[as.character(bd.tab$name),], .) %>% write.csv(., paste("../../Results/Shavings.HubOTUS.", cur.s, ".csv", sep=""), quote=F, row.names=T)
  
  # matches w top 40 RF OTUs
  #bd.tab <- read.csv(paste("../../Results/Shavings.HubOTUS.", cur.s, ".csv", sep=""))
  otunames.rf <- read.csv("../../Results/RF.shavings.t.table.csv", row.names=1)
  x<- intersect(as.character(bd.tab$name), rownames(otunames.rf))
  otunames.rf[x,] %>% write.csv(., paste("../../Results/Shavings.HubOTUS.RF.consensus.", cur.s, ".csv", sep=""), quote=F, row.names=T)

  setwd("../..")

}

## calculate network complexity for each network at each level of p and r

out <-NULL

for (p in c("0.05","0.01","0.005", "0.001")) {
  for (r in c("0.6","0.8")) {
  	out <- rbind(out, cbind(data.frame(State = "IN", P=p, R=r), analyze.network(p, r, "IN", subs=TRUE)))
  	out <- rbind(out, cbind(data.frame(State = "TN", P=p, R=r), analyze.network(p, r, "TN", subs=TRUE)))
  	out <- rbind(out, cbind(data.frame(State = "WA", P=p, R=r), analyze.network(p, r, "WA", subs=TRUE)))
}}
#out$P <- c(rep("0.05",6),rep("0.01",6),rep("0.005",6),rep("0.001",6))
#out$P <- as.factor(out$P)
out$P <- out$P %>% as.character() %>% as.numeric()

out %>% write.csv("Results/Network.complexity.caulosphere.csv")

save.image("R_Environments/Caulo.net.stats.RData")

quartz()
with(out[out$State=="IN" & out$R == 0.8,], plot(P, Kolmogorov, lty=1, type="b", ylim=c(0,max(out$Kolmogorov))))
with(out[out$State=="TN" & out$R == 0.8,], lines(P, Kolmogorov, lty = 2, type="b"))
with(out[out$State=="WA" & out$R == 0.8,], lines(P, Kolmogorov, lty = 3, type="b"))
with(out[out$State=="IN" & out$R == 0.6,], lines(P, Kolmogorov, lty = 1, type="b", pch=2))
with(out[out$State=="TN" & out$R == 0.6,], lines(P, Kolmogorov, lty = 2, type="b", pch=2))
with(out[out$State=="WA" & out$R == 0.6,], lines(P, Kolmogorov, lty = 3, type="b", pch=2))
legend('topleft', legend=c("IN","TN","WA","R=0.8","R=0.6"), lty=c(1,2,3,NA,NA), pch=c(NA,NA,NA,1,2))

quartz()
with(out[out$State=="IN" & out$R == 0.8,], plot(P, Entropy, lty=1, type="b", ylim=c(min(out$Entropy),max(out$Entropy))))
with(out[out$State=="TN" & out$R == 0.8,], lines(P, Entropy, lty = 2, type="b"))
with(out[out$State=="WA" & out$R == 0.8,], lines(P, Entropy, lty = 3, type="b"))
with(out[out$State=="IN" & out$R == 0.6,], lines(P, Entropy, lty = 1, type="b", pch=2))
with(out[out$State=="TN" & out$R == 0.6,], lines(P, Entropy, lty = 2, type="b", pch=2))
with(out[out$State=="WA" & out$R == 0.6,], lines(P, Entropy, lty = 3, type="b", pch=2))
legend('bottomright', legend=c("IN","TN","WA","R=0.8","R=0.6"), lty=c(1,2,3,NA,NA), pch=c(NA,NA,NA,1,2))

quartz()
with(out[out$State=="IN" & out$R == 0.8,], plot(P, Complex.Index, lty=1, type="b", ylim=c(min(out$Complex.Index),max(out$Complex.Index))))
with(out[out$State=="TN" & out$R == 0.8,], lines(P, Complex.Index, lty = 2, type="b"))
with(out[out$State=="WA" & out$R == 0.8,], lines(P, Complex.Index, lty = 3, type="b"))
with(out[out$State=="IN" & out$R == 0.6,], lines(P, Complex.Index, lty = 1, type="b", pch=2))
with(out[out$State=="TN" & out$R == 0.6,], lines(P, Complex.Index, lty = 2, type="b", pch=2))
with(out[out$State=="WA" & out$R == 0.6,], lines(P, Complex.Index, lty = 3, type="b", pch=2))
legend('topleft', legend=c("IN","TN","WA","R=0.8","R=0.6"), lty=c(1,2,3,NA,NA), pch=c(NA,NA,NA,1,2))

quartz()
with(out[out$State=="IN" & out$R == 0.8,], plot(P, Norm.Edge.Complex, lty=1, type="b", ylim=c(min(out$Norm.Edge.Complex),max(out$Norm.Edge.Complex))))
with(out[out$State=="TN" & out$R == 0.8,], lines(P, Norm.Edge.Complex, lty = 2, type="b"))
with(out[out$State=="WA" & out$R == 0.8,], lines(P, Norm.Edge.Complex, lty = 3, type="b"))
with(out[out$State=="IN" & out$R == 0.6,], lines(P, Norm.Edge.Complex, lty = 1, type="b", pch=2))
with(out[out$State=="TN" & out$R == 0.6,], lines(P, Norm.Edge.Complex, lty = 2, type="b", pch=2))
with(out[out$State=="WA" & out$R == 0.6,], lines(P, Norm.Edge.Complex, lty = 3, type="b", pch=2))
legend('topleft', legend=c("IN","TN","WA","R=0.8","R=0.6"), lty=c(1,2,3,NA,NA), pch=c(NA,NA,NA,1,2))

quartz()
with(out[out$State=="IN" & out$R == 0.8,], plot(P, Shannon.ent, lty=1, type="b", ylim=c(min(out$Shannon.ent),max(out$Shannon.ent))))
with(out[out$State=="TN" & out$R == 0.8,], lines(P, Shannon.ent, lty = 2, type="b"))
with(out[out$State=="WA" & out$R == 0.8,], lines(P, Shannon.ent, lty = 3, type="b"))
with(out[out$State=="IN" & out$R == 0.6,], lines(P, Shannon.ent, lty = 1, type="b", pch=2))
with(out[out$State=="TN" & out$R == 0.6,], lines(P, Shannon.ent, lty = 2, type="b", pch=2))
with(out[out$State=="WA" & out$R == 0.6,], lines(P, Shannon.ent, lty = 3, type="b", pch=2))
legend('topleft', legend=c("IN","TN","WA","R=0.8","R=0.6"), lty=c(1,2,3,NA,NA), pch=c(NA,NA,NA,1,2))

quartz()
with(out[out$State=="IN" & out$R == 0.8,], plot(P, Shannon.2nd, lty=1, type="b", ylim=c(min(out$Shannon.2nd),max(out$Shannon.2nd))))
with(out[out$State=="TN" & out$R == 0.8,], lines(P, Shannon.2nd, lty = 2, type="b"))
with(out[out$State=="WA" & out$R == 0.8,], lines(P, Shannon.2nd, lty = 3, type="b"))
with(out[out$State=="IN" & out$R == 0.6,], lines(P, Shannon.2nd, lty = 1, type="b", pch=2))
with(out[out$State=="TN" & out$R == 0.6,], lines(P, Shannon.2nd, lty = 2, type="b", pch=2))
with(out[out$State=="WA" & out$R == 0.6,], lines(P, Shannon.2nd, lty = 3, type="b", pch=2))
legend('topleft', legend=c("IN","TN","WA","R=0.8","R=0.6"), lty=c(1,2,3,NA,NA), pch=c(NA,NA,NA,1,2))

## hubs shared between states

in.hubs <- read.csv("Results/Shavings.HubOTUS.IN.csv")
tn.hubs <- read.csv("Results/Shavings.HubOTUS.TN.csv")
wa.hubs <- read.csv("Results/Shavings.HubOTUS.WA.csv")

intersect(in.hubs$name, tn.hubs$name)
intersect(wa.hubs$name, tn.hubs$name)
intersect(wa.hubs$name, in.hubs$name)
intersect(wa.hubs$name, in.hubs$name) %>% intersect(tn.hubs$name)

library(VennDiagram)

venn.diagram(x = list(in.hubs$name, wa.hubs$name, tn.hubs$name), category.names=c("IN","WA","TN"), filename="Figures/Shavings.Net.Venn.png", output=T)
 
## figures

figure.network <- function (p, r, s, subs=F) {

  gr <- read.csv(paste("Net.analysis/CorrNetworks.", s,"/cyto.shavings.merged.p.", p, ".r.", r,".txt", sep = ''), sep='\t')[,1:3] %>% as.matrix() %>% graph_from_data_frame(directed =F)
  
  edge_attr(gr, name= 'r.sign') <- read.csv(paste("Net.analysis/CorrNetworks.", s,"/cyto.shavings.merged.p.", p, ".r.", r,".txt", sep = ''), sep='\t')[,3] %>% sign() %>% as.factor()
  edge_attr(gr, name= 'r.abs') <- read.csv(paste("Net.analysis/CorrNetworks.", s,"/cyto.shavings.merged.p.", p, ".r.", r,".txt", sep = ''), sep='\t')[,3] %>% abs() %>% as.numeric()
  vertex_attr(gr, name='deg') <- degree(gr)
  vertex_attr(gr, name='fu.ba') <- (vertex_attr(gr, 'name') %>% grepl('Botu',.) ) %>% as.numeric()
  
  bd.tab <- read.csv(paste("Results/Shavings.HubOTUS.", s, ".csv", sep=""), row.names = 1)
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
  
  layout <- layout_with_kk(gr, weights=2-as.numeric(edge_attr(gr, 'r')))
  
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

for (p in c("0.05","0.01","0.005", "0.001")) {
  for (r in c("0.6","0.8")) {
    
    pdf(paste("Net.analysis/Net.figures/Shavings/Allstates.networks.caulo.", "p", p, ".r", r, ".nosubnets.pdf", sep=""), width=12, height=6)
    
    grid.arrange(
      figure.network(p, r, "IN"),
      figure.network(p, r, "TN"),
      figure.network(p, r, "WA"),
      nrow=1,
      top=paste("p = ", p, "r = ", r, sep="")
    )
    
    dev.off()
    
    pdf(paste("Net.analysis/Net.figures/Shavings/Subnetworks/Allstates.networks.caulo.", "p", p, ".r", r, ".wsubnets.pdf", sep=""), width=12, height=6)
    
    grid.arrange(
      figure.network(p, r, "IN", T),
      figure.network(p, r, "TN", T),
      figure.network(p, r, "WA", T),
      nrow=1,
      top=paste("p = ", p, "r = ", r, sep="")
    )
    
    dev.off()

  }}