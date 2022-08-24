
### Network Objects
library(reshape2)
library(sna)
library(network)
library(tnet)
library(igraph)
library(statnet)
#library(bipartite)
library(dplyr)
#library(GGally)
#library(intergraph)
#library(ggnetwork)
library(tibble)
#library(brainGraph)
#library(ggbipart)



# Two-mode Networks
edge <- read.csv("edge_list.csv")

e <- edge[c(1,2)]

# I am building adjacency matrices here so that I can read the weights more easily
tb <- as.data.frame(table(e))
bi.adj0 <- dcast(tb, QID ~ lang, value.var = "Freq")
bi.adj0[is.na(bi.adj0)] <- 0
bi.adj0[,order(colnames(bi.adj0))]
bi.adj <- bi.adj0[,-1]
rownames(bi.adj) <- bi.adj0[,1]


## Projection, p2p
p2p.adj=as.matrix(bi.adj) %*% t(bi.adj)
diag(p2p.adj) <- 0  ## removing self-loops
plot.sociomatrix(p2p.adj,cex.lab = 0.2,font.lab=2)

mean(p2p.adj)
median(p2p.adj)
max(p2p.adj)

hist(c(p2p.adj),ylab="Frequency",xlab="No. of Wikis",main="Degree Distribution of Shared Wikis")


#One-mode projection:wikis
w2w.adj=t(bi.adj) %*% as.matrix(bi.adj) 
diag(w2w.adj) <- 0  ## removing self-loops
plot.sociomatrix(w2w.adj,cex.lab = 0.3,font.lab=2)

mean(w2w.adj)
median(w2w.adj)
max(w2w.adj)

hist(c(w2w.adj),ylab="Frequency",xlab="No. of Shared Policies",main="Degree Distribution of Shared Policies")


#Projection in Igraph  --- not showing degrees/weights
#proj_g <- bipartite.projection(bi.g)


#plot(proj_g$proj1,main="One-mode Projection: Policies",
#     vertex.size=5, vertex.label.cex=0.5,edge.curved=1)
#plot(proj_g$proj2,main="One-mode Projection: Wikis",
#     vertex.size=4, vertex.label.cex=0.3,edge.curved=1)

### descriptive statistics
network.density(bi.net)  #density 0.07791199
graph.motifs.no(g, size = 4) #9793441 cycles

## centrality
# degree centrality for p2p
p2p.degree <- degree_w(net=p2p.adj, measure=c("degree","output","alpha"), alpha=0.5)
hist(as.data.frame(p2p.degree)$alpha)
# closeness centrality for p2p
p2p.closeness <- closeness_w(p2p.adj, alpha=0.5)
hist(as.data.frame(p2p.closeness)$closeness)

#Igraph Object for p2p
p2p.g=graph_from_adjacency_matrix(p2p.adj, mode="undirected", weighted=T)
E(p2p.g)$weight

plot(p2p.g, edge.width=E(p2p.g)$weight/4,
     vertex.size=5, vertex.label.cex=0.5,edge.curved=1,
     layout=layout_randomly(p2p.g))

p2p.eigen <- page_rank(p2p.g, directed = F, damping = 0.85, weights = edge_attr(p2p.g, "weight"))
hist(p2p.eigen$vector)

# merging centrality measures
p2p.ctr <- tibble::rownames_to_column(as.data.frame(p2p.eigen$vector), "QID")
p2p.centrality <- cbind(p2p.degree,p2p.closeness,p2p.ctr)
#write.csv(p2p.centrality, "p2p.centrality.csv")

# degree centrality for w2w
w2w.degree <- degree_w(net=w2w.adj, measure=c("degree","output","alpha"), alpha=0.5)
hist(as.data.frame(w2w.degree)$alpha)
# closeness centrality for w2w
w2w.closeness <- closeness_w(w2w.adj, alpha=0.5)
hist(as.data.frame(w2w.closeness)$closeness)
## plot
w2w.g=graph_from_adjacency_matrix(w2w.adj, mode="undirected", weighted=T)
plot(w2w.g, edge.width=E(w2w.g)$weight/2,
     vertex.size=5, vertex.label.cex=0.5,edge.curved=1,
     vertex.label.family="Arial Bold",
     vertex.color="orange", vertex.frame.color="orange",
     layout=layout_with_gem(w2w.g))+
  geom_nodetext_repel()


# eigenvector centrality for w2w
w2w.eigen <- page_rank(w2w.g, directed = F, damping = 0.85, weights = edge_attr(w2w.g, "weight"))
hist(w2w.eigen$vector)

# merging centrality measures for w2w
w2w.ctr <- tibble::rownames_to_column(as.data.frame(w2w.eigen$vector), "Wiki")
w2w.centrality <- cbind(w2w.degree,w2w.closeness,w2w.ctr)
#write.csv(w2w.centrality, "w2w.centrality.csv")

## edge list for Gephi
edge.p2p = cbind(get.edgelist(p2p.g), E(p2p.g)$weight)
colnames(edge.p2p) = c("Source","Target", "Weight")
#write.csv(edge.p2p,"p2p.edgelist.csv")
#edgelist_policy2policy.csv
edge.w2w = cbind(get.edgelist(w2w.g), E(w2w.g)$weight)
colnames(edge.w2w) = c("Source","Target","Weight")
#write.csv(edge.w2w,"w2w.edgelist.csv")

#Visualizaing the two-mode

colnames(edge) <- c("node","node")
n <- as.data.frame(unique(rbind(edge[c(1)],edge[c(2)])))
n$type <- ifelse(startsWith(n$node, "Q") == TRUE, "policy","wiki")

bi.net <- as.network(e, 
                     directed = FALSE, 
                     vertices = n, 
                     bipartite = TRUE)

gplot(bi.net, 
      gmode="twomode", 
      usearrows = FALSE, 
      displaylabels = F, 
      edge.col="gray")

get.vertex.attribute(bi.net)
network.vertex.names(bi.net)


#Try igraph visualization
g <- graph_from_edgelist(as.matrix(e), directed=F)

is_bipartite(g)
# [1] FALSE
V(g)$type <- V(g)$name %in% e$QID
is_bipartite(g)
V(g)$type

plot(g,main="Policy-Wiki Two-mode Network",
     vertex.size=5, vertex.label.cex=0.5,edge.curved=1)

#detach("package:igraph", unload=TRUE)


# Block models
load("block.fit.Rdata")
# Block models for p2p graph
#fit two-class blockmodels #
ct=list()    # Core w/in,out ties
c=list()      # Isolated core
bic=list()     # Two cores
no_c=list()     # Bipartite structure - no cores!
no_cp=list()  # Ignore core/periphery relations

# fitting 

ct=block.fit(p2p.adj,c(1,1,1,0))    # Core w/in,out ties
c=block.fit(p2p.adj,c(1,0,0,0))         # Isolated core
bic=block.fit(p2p.adj,c(1,0,0,1))        # Two cores
no_c=block.fit(p2p.adj,c(0,1,1,0))        # Bipartite structure - no cores!
no_cp=block.fit(p2p.adj,c(1,NA,NA,0))   # Ignore core/periphery relations

# Core w/in,out ties
lab=ct$block.membership[ct$order.vector]
plot.sociomatrix(ct$blocked.data,labels=list(lab,lab), 
                 main = "Core with in and out Ties")

# Isolated core
lab=c$block.membership[c$order.vector]
plot.sociomatrix(c$blocked.data,labels=list(lab,lab),
                 main = "Isolated Core")

# Two cores
lab=bic$block.membership[bic$order.vector]
plot.sociomatrix(bic$blocked.data,labels=list(lab,lab),
                 main = "Two Cores")

# Bipartite structure - no cores!
lab<-no_c$block.membership[no_c$order.vector]
plot.sociomatrix(no_c$blocked.data,labels=list(lab,lab),
                 main = "Bipartite structure")

# Ignore core/periphery relations
lab<-no_cp$block.membership[no_cp$order.vector]
plot.sociomatrix(no_cp$blocked.data,labels=list(lab,lab),
                 main = "Ignore Core-Peripheral Relations")


# 2(c) Examine the goodness of fit
gof_ct=ct$block.gof
gof_c=c$block.gof
gof_bic=bic$block.gof
gof_no_c=no_c$block.gof


gof=rbind(gof_ct,gof_c,gof_bic,gof_no_c)
rownames(gof)=c("core.ties","isolated.core","two.cores","no.core" )
gof

# R-square
gof_sq=gof^2
gof_sq


## isolated core structure fits the best
View(c)
# Isolated core
lab=c$block.membership[c$order.vector]
plot.sociomatrix(c$blocked.data,labels=list(lab,lab),
                 main = "Isolated Core")

View(c$block.membership)
rownames(p2p.adj)
blocks=as.data.frame(cbind(rownames(p2p.adj),as.data.frame(c$block.membership)))
colnames(blocks) = c("QID","block.membership")
core.p2p = merge(blocks, unique(edge[c(1,3)]), by="id", all=F)
write.csv(core.p2p, "core.p2p.csv")

# Block models for w2w graph
#fit two-class blockmodels #
ct=list()    # Core w/in,out ties
c=list()      # Isolated core
bic=list()     # Two cores
no_c=list()     # Bipartite structure - no cores!
no_cp=list()  # Ignore core/periphery relations

# fitting 

ct=block.fit(w2w.adj,c(1,1,1,0))    # Core w/in,out ties
c=block.fit(w2w.adj,c(1,0,0,0))         # Isolated core
bic=block.fit(w2w.adj,c(1,0,0,1))        # Two cores
no_c=block.fit(w2w.adj,c(0,1,1,0))        # Bipartite structure - no cores!
no_cp=block.fit(w2w.adj,c(1,NA,NA,0))   # Ignore core/periphery relations

# Core w/in,out ties
lab=ct$block.membership[ct$order.vector]
plot.sociomatrix(ct$blocked.data,labels=list(lab,lab), 
                 main = "Core with in and out Ties")

# Isolated core
lab=c$block.membership[c$order.vector]
plot.sociomatrix(c$blocked.data,labels=list(lab,lab),
                 main = "Isolated Core")

# Two cores
lab=bic$block.membership[bic$order.vector]
plot.sociomatrix(bic$blocked.data,labels=list(lab,lab),
                 main = "Two Cores")

# Bipartite structure - no cores!
lab<-no_c$block.membership[no_c$order.vector]
plot.sociomatrix(no_c$blocked.data,labels=list(lab,lab),
                 main = "Bipartite structure")

# Ignore core/periphery relations
lab<-no_cp$block.membership[no_cp$order.vector]
plot.sociomatrix(no_cp$blocked.data,labels=list(lab,lab),
                 main = "Ignore Core-Peripheral Relations")


# 2(c) Examine the goodness of fit
gof_ct=ct$block.gof
gof_c=c$block.gof
gof_bic=bic$block.gof
gof_no_c=no_c$block.gof


gof=rbind(gof_ct,gof_c,gof_bic,gof_no_c)
rownames(gof)=c("core.ties","isolated.core","two.cores","no.core" )
gof

# R-square
gof_sq=gof^2
gof_sq


## isolated core structure fits the best
View(c)
# Isolated core
lab=c$block.membership[c$order.vector]
plot.sociomatrix(c$blocked.data,labels=list(lab,lab),
                 main = "Isolated Core")

View(c$block.membership)
rownames(w2w.adj)
blocks=as.data.frame(cbind(rownames(w2w.adj),as.data.frame(c$block.membership)))
colnames(blocks) = c("Wiki","block.membership")

wiki.meta=read.csv("wiki_level_metadata-20220707.csv",stringsAsFactors = F)
core.w2w = merge(blocks, wiki.meta, by="Wiki", all = F)

write.csv(core.w2w,"core.w2w.csv")


#k-core p2p
p2p.bi.adj=ifelse(p2p.adj>34, 1, 0)
p2p.bi = graph_from_adjacency_matrix(p2p.bi.adj, mode="undirected")

#plot(p2p.bi,
#     vertex.size=5, vertex.label.cex=0.5,edge.curved=1,
#     layout=layout_with_fr)

coreness <- graph.coreness(p2p.bi) 
maxCoreness <- max(coreness)

k.core.p  <- cbind(rownames(p2p.adj),as.data.frame(coreness(p2p.bi)))

verticesHavingMaxCoreness <- which(coreness == maxCoreness) 
kcore <- induced.subgraph(p2p.bi,vids=verticesHavingMaxCoreness)

plot(kcore, 
     vertex.label=get.vertex.attribute(kcore,name='vert.names',index=V(kcore)))

k.core.vertices <- rownames(subset(k.core.p, k.core.p[2]==32))

p2p.g2 <- delete_vertices(p2p.g, as.vector(k.core.vertices))

V(p2p.g2)
plot(p2p.g2,edge.width=E(p2p.g2)$weight/4,
     vertex.size=5, vertex.label.cex=0.5,edge.curved=1,)


#k-core w2w
w2w.bi.adj=ifelse(w2w.adj>8, 1, 0)
w2w.bi = graph_from_adjacency_matrix(w2w.bi.adj, mode="undirected")

plot(w2w.bi,
     vertex.size=5, vertex.label.cex=0.5,edge.curved=1,
     layout=layout_with_fr)

coreness <- graph.coreness(w2w.bi) 
maxCoreness <- max(coreness)

k.core.w  <- cbind(rownames(w2w.adj),as.data.frame(coreness(w2w.bi)))

verticesHavingMaxCoreness <- which(coreness == maxCoreness) 
kcore <- induced.subgraph(w2w.bi,vids=verticesHavingMaxCoreness)

plot(kcore, 
     vertex.label=get.vertex.attribute(kcore,name='vert.names',index=V(kcore)))

k.core.vertices <- rownames(subset(k.core.w, k.core.w[2]==76))

w2w.g2 <- delete_vertices(w2w.g, as.vector(k.core.vertices))
w2w.bi2 <- delete_vertices(w2w.bi, as.vector(k.core.vertices))

V(w2w.g2)
plot(w2w.g2,edge.width=E(w2w.g2)$weight/4,
     vertex.size=5, vertex.label.cex=0.5,edge.curved=1,)

V(w2w.bi2)
plot(w2w.bi)



#nds_nl ego network
nds.nl=as.data.frame(w2w.adj)
nds.nl = as.data.frame(cbind(rownames(w2w.adj),ndl.nl$`nds-nl`))
nds.nl.ego1 = subset(nds.nl,nds.nl$V2 > 0)
nds.nl.ego2 = subset(nds.nl,nds.nl$V2 > 1)

#hatian ego network
ht=as.data.frame(w2w.adj)
ht = as.data.frame(cbind(rownames(w2w.adj),ndl.nl$`ht`))
ht.ego = subset(ht,ht$V2 > 1)

#keep = append(as.vector(ht.ego$V1),"ht")

sum <- apply(w2w.adj[,c(1:245)], 1, sum)

ego.extract(w2w.net)

#gplot.layout.adj(ego.extract(w2w.net)$`zh-yue`,layout.par = NULL)
gplot(ego.extract(w2w.net)$`jam`,usearrows = FALSE, displaylabels = T,
      vertex.col = "black",
      edge.col = "grey")

jam.ego=graph_from_adjacency_matrix(ego.extract(w2w.net)$`jam`, mode="undirected", weighted=T)
plot(jam.ego, edge.width=E(w2w.g)$weight/3,
     vertex.size=5, vertex.label.cex=1,edge.curved=1,
     vertex.label.family="Arial Black",
     vertex.color="orange", vertex.frame.color="#ffffff",
     layout=layout_with_gem(jam.ego))




#Connectivity of weighted networks
#p2p.bi.adj=ifelse(p2p.adj>1, 1, 0)
p2p.bi.adj=ifelse(p2p.adj>1, 1, 0)
p2p.bi = graph_from_adjacency_matrix(p2p.bi.adj, mode="undirected")
plot(p2p.bi)
vertex_connectivity(p2p.bi, source = NULL, target = NULL, checks = TRUE)


Comp = clusters(p2p.bi)
Comp$csize[Comp$membership]


## S3 method for class 'igraph'
cohesion(p2p.bi, checks = TRUE)


#k-core p2p weighted
s_core(p2p.g, W = "weight")
s.core <- cbind(rownames(p2p.adj),as.data.frame(s_core(p2p.g)))


