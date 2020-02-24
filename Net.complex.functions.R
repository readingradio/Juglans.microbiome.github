library(igraph)
library(tidygraph)
library(graph)
library(QuACN)
library(acss)

setwd("OACC")
source("scripts/BDM2D.R")
setwd("..")

analyze.network <- function (p, r, s, subs=F, bs=4, tissue="shavings") {

  gr <- read.csv(paste("Net.analysis/CorrNetworks.", s,"/cyto.", tissue, ".merged.p.", p, ".r.", r,".txt", sep = ''), sep='\t')[,1:3] %>% as.matrix() %>% graph_from_data_frame(directed =F)
  
  # calculate Shannon entropy and Kolmogorov-Chaitin algorithmic complexity
  
  # For K-C AC, get the adjacency matrix from graph
  adj.mat <- as_adj(gr) %>% as.matrix()
 
  # calculate K
  dims<-(dim(adj.mat) %/% bs) * bs
  K<-bdm2D(adj.mat[1:dims[1],1:dims[2]], blockSize=bs)

  g<-as(as_adj(gr), "graphNEL")
  E<-topologicalInfoContent(g)$entropy
  C<-graphIndexComplexity(g)
  N<-normalizedEdgeComplexity(g)

  # decompose into string and calculate Shannon Entropy
  adj.str <-  adj.mat %>% as.matrix() %>% as.vector() %>% paste(collapse="")
  S <- entropy(adj.str)[[1]]
  S2<- entropy2(adj.str)[[1]]
  
  data.frame(Kolmogorov = K, Entropy = E, Complex.Index = C, Norm.Edge.Complex = N, Shannon.ent = S, Shannon.2nd = S2)
}