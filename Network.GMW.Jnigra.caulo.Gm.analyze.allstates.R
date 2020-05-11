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
library(phyloseq)
library(vegan)
library(reshape2)
library(tanagR)

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
  	
  	IN <- read.network (p, r, "IN", subs=T, tissue="shavings", path = ".")
	TN <- read.network (p, r, "TN", subs=T, tissue="shavings", path = ".")
  	WA <- read.network (p, r, "WA", subs=T, tissue="shavings", path = ".")
  	
  	states <- list(IN,TN,WA)
  	names(states) <- c("IN", "TN", "WA")
  	out <- rbind(out, cbind(data.frame(P = rep(p, 3), R = rep(r, 3)), compare.networks(states)))
}}

out$P <- out$P %>% as.character() %>% as.numeric()

out %>% write.csv("Results/Network.complexity.caulosphere.csv")

#out <- read.csv("Results/Network.complexity.caulosphere.csv")[,-1]
save.image("R_Environments/Caulo.net.stats.RData")

#load("R_Environments/Caulo.net.stats.RData")

quartz()
net.complex.multipanel (out)
net.complex.multipanel (out[out$P<0.05,])

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

# function to create a network object

make.network.graph <- function (p, r, s, subs=F) {

  gr <- read.csv(paste("Net.analysis/CorrNetworks.", s,"/cyto.shavings.merged.p.", p, ".r.", r,".txt", sep = ''), sep='\t')[,1:3] %>% as.matrix() %>% graph_from_data_frame(directed =F)
  
  edge_attr(gr, name= 'r.sign') <- read.csv(paste("Net.analysis/CorrNetworks.", s,"/cyto.shavings.merged.p.", p, ".r.", r,".txt", sep = ''), sep='\t')[,3] %>% sign() %>% as.factor()
  edge_attr(gr, name= 'r.abs') <- read.csv(paste("Net.analysis/CorrNetworks.", s,"/cyto.shavings.merged.p.", p, ".r.", r,".txt", sep = ''), sep='\t')[,3] %>% abs() %>% as.numeric()
  vertex_attr(gr, name='deg') <- igraph::degree(gr)
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
  
  gr
}

## create object for chosen caulosphere graph (can switch between TN and WA)

top.clusters <- function (p='0.01', r='0.8', state) {
	g.r0.8.p0.01 <- make.network.graph(p, r, state)

	# remove negative connections
	# as.factor(sign(c(1,-8)))

	negs<-edge_attr(g.r0.8.p0.01, name= 'r.sign') == 1
	gr.r8p01.pos <-delete_edges(g.r0.8.p0.01, which(negs))

	# pull out major subnetworks

	dec.subrs.r8p01 <- decompose.graph(gr.r8p01.pos)
	clusters(gr.r8p01.pos)
	twobiggest <- rev(order(clusters(gr.r8p01.pos)$csize))[1:2]
	cluster.1 <- dec.subrs.r8p01[[twobiggest[1]]]
	cluster.2 <- dec.subrs.r8p01[[twobiggest[2]]]
	
	vertex_attr(g.r0.8.p0.01, "Cluster") <- (0 + 1*(clusters(gr.r8p01.pos)$membership == twobiggest[1]) + 2*(clusters(gr.r8p01.pos)$membership == twobiggest[2])) %>% as.factor()

	edge_attr(g.r0.8.p0.01, "ClusterEdge") <- (0 + 1*(attr(E(g.r0.8.p0.01),"vnames") %in% attr(E(cluster.1),"vnames")) + 2*(attr(E(g.r0.8.p0.01),"vnames") %in% attr(E(cluster.2),"vnames"))) %>% as.factor()
	
	list(cluster.1=cluster.1, cluster.2=cluster.2, net.annot=g.r0.8.p0.01)

}

wa12 <- top.clusters(state='WA')
tn12 <- top.clusters(state='TN')

# just hub OTUs from the two clusters 
wa.clust.1.otus <- with(wa12,vertex_attr(cluster.1,"name")[as.logical(vertex_attr(cluster.1,"hub"))])
wa.clust.2.otus <- with(wa12,vertex_attr(cluster.2,"name")[as.logical(vertex_attr(cluster.2,"hub"))])

tn.clust.1.otus <- with(tn12,vertex_attr(cluster.1,"name")[as.logical(vertex_attr(cluster.1,"hub"))])
tn.clust.2.otus <- with(tn12,vertex_attr(cluster.2,"name")[as.logical(vertex_attr(cluster.2,"hub"))])

# Trees that tested positive for GM with the probe were WA55RN1, WA130RN4, and WA132RN9. There were five other trees that only tested positive once. Those include the WA132 (the remaining two), WA272RN10, and WAWTBNL22, WAWTBNL23
# rownames(shavings.merge.df[,clust.1.otus])

wa.gm <- c("WA_RN1_55s_S44","WA_RN4_130s_S21","WA_RN9_132s_S61","WA_RN8_132s_S53","WA_RN7_132s_S37","WA_BNL23_WTs_S6","WA_BNL22_WTs_S45")
wa.neg<- c("WA_BNL18_272s_S85","WA_BNL19_55s_S13","WA_BNL21_WTs_S93")

# Trees in TCD positive and negative location in TN

tn.tcd <- rownames(shavings.merge.df)[regexpr("TN_LS",rownames(shavings.merge.df),perl=T)==1]
tn.neg <- rownames(shavings.merge.df)[regexpr("TN",rownames(shavings.merge.df),perl=T)==1] %>% setdiff(tn.tcd)

indian <- rownames(shavings.merge.df)[regexpr("IN_",rownames(shavings.merge.df),perl=T)==1] %>% setdiff(tn.tcd)

intersect(wa.clust.1.otus,tn.clust.1.otus)
intersect(wa.clust.1.otus,tn.clust.2.otus)
intersect(wa.clust.2.otus,tn.clust.1.otus)
intersect(wa.clust.2.otus,tn.clust.2.otus)

shavings.merge.df[] <- lapply(shavings.merge.df, function(x) as.numeric(as.character(x)))

# cluster 1 and 2 for gm
shavings.merge.df[wa.gm,wa.clust.1.otus] %>% t() %>% cbind(tax_table(shavings.merge)[wa.clust.1.otus],.) %>% as.data.frame()
shavings.merge.df[wa.gm,wa.clust.2.otus] %>% t() %>% cbind(tax_table(shavings.merge)[wa.clust.2.otus],.) %>% as.data.frame()

# cluster 1 and 2 for neg trees
shavings.merge.df[wa.neg,wa.clust.1.otus] %>% t() %>% cbind(tax_table(shavings.merge)[wa.clust.1.otus],.) %>% as.data.frame()
shavings.merge.df[wa.neg,wa.clust.2.otus] %>% t() %>% cbind(tax_table(shavings.merge)[wa.clust.2.otus],.) %>% as.data.frame()

# cluster 1  for tn knox and tn polk
shavings.merge.df[tn.tcd,tn.clust.1.otus] %>% t() %>% cbind(tax_table(shavings.merge)[tn.clust.1.otus],.) %>% as.data.frame()
shavings.merge.df[tn.neg,tn.clust.1.otus] %>% t() %>% cbind(tax_table(shavings.merge)[tn.clust.1.otus],.) %>% as.data.frame()

# cluster 2 for tn tn knox and tn polk
shavings.merge.df[tn.tcd,tn.clust.2.otus] %>% t() %>% cbind(tax_table(shavings.merge)[tn.clust.2.otus],.) %>% as.data.frame()
shavings.merge.df[tn.neg,tn.clust.2.otus] %>% t() %>% cbind(tax_table(shavings.merge)[tn.clust.2.otus],.) %>% as.data.frame()

# combine wa cluster 1 and 2, tn cluster 1 and 2

hm <- shavings.merge.df[c(indian,tn.neg,tn.tcd,wa.neg,wa.gm),c(tn.clust.1.otus,tn.clust.2.otus,wa.clust.1.otus,wa.clust.2.otus)] %>% t() %>% cbind(tax_table(shavings.merge)[c(tn.clust.1.otus,tn.clust.2.otus,wa.clust.1.otus,wa.clust.2.otus)],.) %>% as.data.frame()

hm[,-c(1:7)] <- lapply(hm[,-c(1:7)], function(x) as.numeric(as.character(x)))

#hm.melt <- (1-(hm[,-c(1:7)] %>% decostand(method="log") %>% decostand(.,method="max", MARGIN=1))) %>% cbind(data.frame(OTU=rownames(hm), Cluster = c(rep(1,length(tn.clust.1.otus)),rep(2,length(tn.clust.2.otus)),rep(3,length(wa.clust.1.otus)),rep(4,length(wa.clust.2.otus))))) %>% melt(id.vars=c("OTU","Cluster"),variable.name="Sample") 

hm.melt <- hm[,-c(1:7)] %>% decostand(method="log") %>% decostand(.,method="max", MARGIN=1) %>% cbind(data.frame(OTU=rownames(hm), Cluster = as.factor(c(rep(1,length(tn.clust.1.otus)),rep(2,length(tn.clust.2.otus)),rep(3,length(wa.clust.1.otus)),rep(4,length(wa.clust.2.otus)))))) %>% melt(id.vars=c("OTU","Cluster"),variable.name="Sample") 

hm.melt$Sample <- factor(hm.melt$Sample, c(indian,tn.neg,tn.tcd,wa.neg,wa.gm))

hm.melt$OTU <- factor(hm.melt$OTU, rev(c(tn.clust.1.otus,tn.clust.2.otus,wa.clust.1.otus,wa.clust.2.otus)))

axislabels <- with(hm, paste(rownames(hm), gsub("[sgfocp]__","",Genus,perl=T) %>% gsub("_"," ",.,perl=T) %>% gsub("[\\s]?(unclassified|ge)[\\s]?","",.,perl=T)%>% gsub("\\([A-Za-z0-9 _]*\\)"," ",.,perl=T)))
names(axislabels) <- rownames(hm)

i <- grep("(uncultured$)|([01234567890-]+$)|(group$)", axislabels, perl=T)
axislabels[i]<-paste(names(axislabels)[i], hm$Family[i])
j <- grep("(uncultured$)|([01234567890-]+$)|(group$)|([01234567890-]+_fa$)", axislabels, perl=T)
axislabels[j]<-paste(names(axislabels)[j], hm$Order[j])
k <- grep("(uncultured_fa$)|([01234567890-]+_or$)", axislabels, perl=T)
axislabels[k]<-paste(names(axislabels)[k], hm$Phylum[k])
l <- grep("[01234567890-]+$", axislabels, perl=T)
axislabels[l]<-paste(names(axislabels)[l], hm$Class[l])

samplelabels<-c(rep("Indiana", length(indian)),rep("TN Polk Co. Clones", length(tn.neg)-3),rep("TN Polk Co. Bedrun", 3),rep("TN Knox Co. Wild", length(tn.tcd)),rep("WA TCD-negative", length(wa.neg)),rep("WA TCD-positive", length(wa.gm)))

colors.subnetworks<-tanagr_palette("dacnis_berlepschi")[1:4]

hm.melt %>% ggplot(., aes(x=Sample, y=OTU, fill=value)) + geom_tile() + scale_fill_gradient(low="white", high="black") + theme(axis.text.x = element_text(angle=90, size=6), legend.position = "none", axis.text.y = element_text(hjust=0))+ scale_y_discrete(labels=axislabels)

heatmap.plot <- hm.melt %>% ggplot(., aes(x=Sample, y=OTU, fill=Cluster, alpha=value)) + geom_tile() + scale_fill_manual(values=colors.subnetworks) + scale_alpha_continuous(range=c(0,1)) + scale_y_discrete(labels=axislabels) + scale_x_discrete(labels=samplelabels)+
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, size=6),
  legend.position = "none",
  axis.text.y = element_text(hjust=0),
  axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()) 

## function to create netwrok graph object

figure.network <- function (p, r, s, subs=F, subcolors=NULL, annot.obj=NULL) {

  if (is.null(annot.obj)) gr <- make.network.graph(p, r, s, subs=F, subcolors=subcolors)
  else gr <- annot.obj
  
  layout <- layout_with_kk(gr, weights=2-as.numeric(edge_attr(gr, 'r')))
  
  if(is.null(subcolors))
	ggraph(gr, layout) + 
    geom_edge_link(aes(width=.1, alpha=r.abs, color=as.factor(r.sign))) + 
    scale_edge_color_manual(values=c('grey','blue')) + 
    scale_edge_width(range=c(0.2,1)) + 
    geom_node_point(aes(size=deg, shape=as.factor(hub), color=as.factor(fu.ba))) + 
    scale_size(range = c(.2,2.5)) + scale_color_manual(values=c('red','black','green')) +
    scale_shape_manual(values=c(1,19,8)) + 
    theme_graph() + 
    theme(legend.position='none')
  else
	ggraph(gr, layout) + 
    geom_edge_link(aes(width=.1, alpha=r.abs, color=as.factor(ClusterEdge))) + 
    scale_edge_color_manual(values=c('grey',subcolors)) + 
    scale_edge_width(range=c(0.2,1)) + 
    geom_node_point(aes(size=deg, shape=as.factor(hub), color=as.factor(Cluster))) + 
    scale_size(range = c(.2,2.5)) + scale_color_manual(values=c('grey',subcolors)) +
    scale_shape_manual(values=c(1,19,8)) + 
    theme_graph() + 
    theme(legend.position='none')
}

## make side networks as a legend

<- figure.network(subcolors = colors.subnetworks[1:2], annot.obj=tn12$net.annot)
<- figure.network(subcolors = colors.subnetworks[3:4], annot.obj=wa12$net.annot)
heatmap.plot

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