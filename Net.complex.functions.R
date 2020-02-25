library(tidyverse)
library(igraph)
library(tidygraph)
library(graph)
library(QuACN)
library(acss)
library(doParallel)
library(foreach)
library(doBy)

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

net.complex.multipanel <- function (x) {
	par(mfrow=c(2,3))
	for (i in c("Kolmogorov", "Entropy", "Complex.Index", "Norm.Edge.Complex", "Shannon.ent")) network.complexity.plot(x, i)
	network.complexity.plot(x, "Shannon.2nd", leg="topleft")
}