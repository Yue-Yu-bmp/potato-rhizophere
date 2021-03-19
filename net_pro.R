net_pro<-function(igraph,outdir){

  num.edges <- length(E(igraph)) 
  num.edges

  num.vertices <- length(V(igraph))
  num.vertices

  connectance <- edge_density(igraph,loops=FALSE)
  connectance

  average.degree <- mean(igraph::degree(igraph))
  average.degree
  
  average.path.length <- average.path.length(igraph)
  average.path.length

  diameter <- diameter(igraph, directed = FALSE, unconnected = TRUE, weights = NA)
  diameter

  edge.connectivity <- edge_connectivity(igraph)
  edge.connectivity

  clustering.coefficient <- transitivity(igraph) 
  clustering.coefficient
  no.clusters <- no.clusters(igraph)
  no.clusters

  centralization.degree <- centralization.degree(igraph)$centralization
  centralization.degree

  centralization.betweenness <- centralization.betweenness(igraph)$centralization 
  centralization.betweenness

  centralization.closeness <- centralization.closeness(igraph)$centralization
  centralization.closeness
  
  num.pos.edges<-sum(E(igraph)$weight>0)# number of postive correlation
  num.neg.edges<-sum(E(igraph)$weight<0)# number of negative correlation
  
  igraph.network.pro <- rbind(num.edges,num.pos.edges,num.neg.edges,num.vertices,connectance,average.degree,average.path.length,diameter,edge.connectivity,clustering.coefficient,no.clusters,centralization.degree,centralization.betweenness,centralization.closeness)
  rownames(igraph.network.pro)<-c("num.edges","num.pos.edges","num.neg.edges","num.vertices","connectance","average.degree","average.path.length","diameter","edge.connectivity","clustering.coefficient","no.clusters","centralization.degree","centralization.betweenness","centralization.closeness")
  colnames(igraph.network.pro)<- "value"
  igraph.network.pro
}