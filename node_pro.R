node_pro<-function(igraph,outdir){  

  igraph.degree<-igraph::degree(igraph)

  igraph.cen.degree<-centralization.degree(igraph)$res

  igraph.betweenness<-centralization.betweenness(igraph)$res

  igraph.closeness<-centralization.closeness(igraph)$res
  
  igraph.node.pro <- cbind(igraph.degree,igraph.closeness,igraph.betweenness,igraph.cen.degree)
  colnames(igraph.node.pro)<-c("igraph.degree","igraph.closeness","igraph.betweenness","igraph.cen.degree")
  igraph.node.pro
}