library(tidyverse)
library(igraph)
library(tidygraph)
library(ggraph)
library(graph)
library(QuACN)
library(acss)
library(doParallel)
library(foreach)
library(doBy)
library(Hmisc)

setwd("OACC")
source("scripts/BDM2D.R")
setwd("..")

read.network <- function (p, r, s, subs=F, tissue="shavings", path=".") read.csv(paste(path, "/Net.analysis/CorrNetworks.", s,"/cyto.", tissue, ".merged.p.", p, ".r.", r,".txt", sep = ''), sep='\t')[,1:3] %>% as.matrix() %>% graph_from_data_frame(directed =F)

analyze.network <- function (gr=NULL, adj.mat=NULL, bs=4) {

  # calculate Shannon entropy and Kolmogorov-Chaitin algorithmic complexity
  
  # For K-C AC, get the adjacency matrix from graph
  if (is.null(adj.mat)) {
  	if (is.null(gr)) return ("Error")
  	adj.mat <- as_adj(gr) %>% as.matrix()
  }

  # calculate K
  dims<-(dim(adj.mat) %/% bs) * bs
  K<-bdm2D(adj.mat[1:dims[1],1:dims[2]], blockSize=bs)

  g<-as(adj.mat, "graphNEL")
  E<-topologicalInfoContent(g)$entropy
  C<-graphIndexComplexity(g)
  N<-normalizedEdgeComplexity(g)

  # decompose into string and calculate Shannon Entropy
  adj.str <-  adj.mat %>% as.matrix() %>% as.vector() %>% paste(collapse="")
  S <- entropy(adj.str)[[1]]
  S2<- entropy2(adj.str)[[1]]
  
  data.frame(Kolmogorov = K, Entropy = E, Complex.Index = C, Norm.Edge.Complex = N, Shannon.ent = S, Shannon.2nd = S2)
}

## compare.networks ()
## This function attempts to correct for network size when network complexity measures among networks
## by randomly subsampling adjacency matrices. It finds the dimension of the smallest network and uses
## this dimension to create a bootstrap sample of network complexity measurements.

compare.networks <- function (g = list(), bs=4, subs=T) { # takes a named list of networks
	
	# first check size of each network
	sizes<-NULL; for (i in 1:length(g)) sizes <- c(sizes, (as_adj(g[[i]]) %>% as.matrix() %>% dim())[1])

	# create size that is the largest dimension of smallest adjacency matrix that is divisible by the block size to be used by Kolmogrov
	cutoff.size <- ( min(sizes) %/% bs ) * bs

	# now do a bootstrap subsample of the network by taking a random contiguous block within the network, and calculate standard deviation
	
	# use parallel processing
	set.seed(84985)
	registerDoParallel(8)
	
	# create a bootstrapped dataframe of complexity measures for randomly selected subgraphs
	boot.out <- foreach (Try=1:20, .combine=rbind) %dopar% {
		out <- NULL
		for (i in 1:length(g)) {
			m <- as_adj(g[[i]]) %>% as.matrix()
			v <- sample (sizes[i])							# create a vector to randomly permute the rows and columns together in the same order
			m <- m[v,v] 										# randomly permute the rows and columns
			start <- sample(sizes[i] - cutoff.size, 1)		# randomly choose what index to start random block of size = cutoff.size from the permuted matrix
			indices <- start:(start + cutoff.size - 1)		# figure out what index to end and make vector of indices to use for subsample
			subm <- m[indices,indices]						# take out random block of size = cutoff.size from the permuted matrix
			out <- rbind(out, cbind(data.frame(Name = (names(g))[i]), analyze.network(adj.mat = subm)))		# analyze randomly subsetted matrix
		}
		out
	}

	# now return mean and standard error for each of the metrics for each matrix by its name
	summaryBy( . ~ Name, boot.out, FUN = c(mean, sd) )
}

network.complexity.plot <- function (x, stat=c("Kolmogorov", "Entropy", "Complex.Index", "Norm.Edge.Complex", "Shannon.ent", "Shannon.2nd"), maxx=max(x$P), maxy=max(x$y2), leg=NA) {

	s.m <- paste(stat, "mean", sep=".")
	s.s <- paste(stat, "sd", sep=".")
	x$y1 <- x[,s.m] - x[,s.s]
	x$y2 <- x[,s.m] + x[,s.s]

	x[x$Name=="IN" & x$R == 0.8, c("P", s.m)] %>% plot(lty=1, type="b", xlim=c(0,maxx), ylim=c(0,maxy))
	x[x$Name=="TN" & x$R == 0.8, c("P", s.m)] %>% lines(lty = 2, type="b")
	x[x$Name=="WA" & x$R == 0.8, c("P", s.m)] %>% lines(lty = 3, type="b")

	x[x$Name=="IN" & x$R == 0.6, c("P", s.m)] %>% lines(lty = 1, type="b", pch=2)
	x[x$Name=="TN" & x$R == 0.6, c("P", s.m)] %>% lines(lty = 2, type="b", pch=2)
	x[x$Name=="WA" & x$R == 0.6, c("P", s.m)] %>% lines(lty = 3, type="b", pch=2)
	
	for (a in c("IN", "TN", "WA")) for (b in c(.8, .6)) for (c in c("y1", "y2")) x[x$Name==a & x$R==b,"P"] %>% arrows(., x[x$Name==a & x$R==b,s.m], ., x[x$Name==a & x$R==b,c], angle=90, length=.025)
	
	legend(leg, legend=c("IN","TN","WA","R=0.8","R=0.6"), lty=c(1,2,3,NA,NA), pch=c(NA,NA,NA,1,2))
}

net.complex.multipanel <- function (x,which.leg="last") {
	par(mfrow=c(2,3))
	if(which.leg=="last"){
	for (i in c("Kolmogorov", "Entropy", "Complex.Index", "Norm.Edge.Complex", "Shannon.ent")) network.complexity.plot(x, i)
	network.complexity.plot(x, "Shannon.2nd", leg="topleft")
	} else if (which.leg=="first") {
	network.complexity.plot(x, "Kolmogorov", leg="topleft")
	for (i in c("Entropy", "Complex.Index", "Norm.Edge.Complex", "Shannon.ent", "Shannon.2nd")) network.complexity.plot(x, i)
	}
}


# function to create a network object

make.network.graph <- function (p, r, s, subs=F, tissue="shavings") {

  gr <- read.csv(paste("Net.analysis/CorrNetworks.", s,"/cyto.", tissue, ".merged.p.", p, ".r.", r,".txt", sep = ''), sep='\t')[,1:3] %>% as.matrix() %>% graph_from_data_frame(directed =F)
  
  edge_attr(gr, name= 'r.sign') <- read.csv(paste("Net.analysis/CorrNetworks.", s,"/cyto.", tissue, ".merged.p.", p, ".r.", r,".txt", sep = ''), sep='\t')[,3] %>% sign() %>% as.factor()
  edge_attr(gr, name= 'r.abs') <- read.csv(paste("Net.analysis/CorrNetworks.", s,"/cyto.", tissue, ".merged.p.", p, ".r.", r,".txt", sep = ''), sep='\t')[,3] %>% abs() %>% as.numeric()
  vertex_attr(gr, name='deg') <- igraph::degree(gr)
  vertex_attr(gr, name='fu.ba') <- (vertex_attr(gr, 'name') %>% grepl('Botu',.) ) %>% as.numeric()

  bd.tab <- read.csv(paste("Results/", capitalize(tissue), ".HubOTUS.", s, ".csv", sep=""), row.names = 1)
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

## create object for chosen caulosphere or soil graph (can switch between TN and WA)

top.clusters <- function (p='0.01', r='0.8', state, tissue="shavings") {
	g.r0.8.p0.01 <- make.network.graph(p, r, state, tissue=tissue)

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

figure.network <- function (p, r, s, subs=F, subcolors=NULL, annot.obj=NULL, tissue="shavings") {

  if (is.null(annot.obj)) gr <- make.network.graph(p, r, s, subs=F, tissue=tissue)
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
    theme(legend.position='none', plot.margin=unit(c(0,0,0,0),"cm"))
  else
	ggraph(gr, layout) + 
    geom_edge_link(aes(width=.1, alpha=r.abs, color=as.factor(ClusterEdge))) + 
    scale_edge_color_manual(values=c('grey',subcolors)) + 
    scale_edge_width(range=c(0.2,1)) + 
    geom_node_point(aes(size=deg, shape=as.factor(hub), color=as.factor(Cluster))) + 
    scale_size(range = c(.2,2.5)) + scale_color_manual(values=c('grey',subcolors)) +
    scale_shape_manual(values=c(1,19,8)) + 
    theme_graph() + 
    theme(legend.position='none', plot.margin=unit(c(0,0,0,0),"cm"))
}

