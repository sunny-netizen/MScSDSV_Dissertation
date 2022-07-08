# code to plot the hierarchical trees
# using size colours according to colours in map for each threshold
#in order to urse ggnet we need the following package
#install.packages("GGally")
#library(GGally)
#install.packages("visNetwork")
#install.packages("networkD3")

library(ggplot2)
library(visNetwork)
library(networkD3)
library(igraph)
library(RColorBrewer)

# Let us select the thresholds used to construct the tree in "hierarchical_tree.R"
#v_jumps <- c(50,80,90,100,110,130,150,165,195,205,210,250,400)
v_jumps <- c(25,26,32, 33)

# path for results
path_res <- paste0("/Users/yun/Documents/CASA/MScSDSV_Dissertation/run_spercolation_outputs/membTables/")
path_plots <- paste0("/Users/yun/Documents/CASA/MScSDSV_Dissertation/tree_outputs/")

#file with the tree network nodes named after levels and clusters
file_tree_graph <- paste0(path_plots,"tree_graph.txt")

#read the file as graph
g_tree <- read_graph(file_tree_graph, format = "ncol")
#let us get the parent node
tail(V(g_tree))
v_dim=length(V(g_tree))
tail(as_edgelist(g_tree))
g_parent=as.character(as_edgelist(g_tree)[v_dim-1,][2])

#let us plot the trees removing nodes with low degree
#or with no children
tail(as_edgelist(g_tree))
V(g_tree)$colour="gray"
#plot(g_tree,vertex.size=log(V(g_tree)$size),vertex.label=NA,layout=layout_as_tree(g_tree,root = g_parent))
plot(g_tree,vertex.size=1,vertex.label=NA,layout=layout_as_tree(g_tree,root = g_parent),
     vertex.color=V(g_tree)$colour)
title("Hierarchical tree, London ",cex.main=2)
  
  deg_ver=degree(g_tree,V(g_tree))
  #max(deg_ver)
  #min(deg_ver)
  #for individual degrees use deg_ver[[i]] for vertex in position i in V(g_tree)[i]
  
  #in order to plot the nodes with the same colours as the maps, we need at each level define the colours again
  top10 = colors()[c(553,29,258,654,91,115,456,48,102,40)]
  
  
  #let us get the characteristics of the nodes
  i_level=1
  for(i_tmp in v_jumps)
  {
    i_tmp=v_jumps[i_level]
    #let us get the n. of intersections for the clusters 
    #this file has already been sorted
    file_n_clusts <- paste0("Results/n_clusters_p",i_tmp,".txt")
    data_clusters <- read.table(file_n_clusts, header = T)
    tail(data_clusters)
    dim(data_clusters)
    dat_clusts <- data_clusters[data_clusters$n_points>=50,]
    tail(dat_clusts)
    dim(dat_clusts)
    d_clusts = dim(dat_clusts)[1]
    
    #firs we need to match the name of the nodes in the tree with the clusters 
    if(i_level>=2)
    {
      dat_clusts$id_Tree=paste0("L",i_level+1,"_",dat_clusts$id_cluster)
    }else{
      dat_clusts$id_Tree=dat_clusts$id_cluster
    }
    head(dat_clusts)
    if(d_clusts<=10)
    {
      #assign colours to top 15
      dat_clusts$colour=top10[1:d_clusts]
      head(dat_clusts)
      
      v=match(V(g_tree)$name,as.character(dat_clusts$id_Tree))
      v_pos_nn=which(!is.na(v))
      V(g_tree)[v_pos_nn]$name
      for(i_pos in 1:length(v_pos_nn))
      {
        #we need to find 
        #i_pos=1
        tmp_data=dat_clusts[dat_clusts$id_Tree==V(g_tree)[v_pos_nn[i_pos]]$name,]
        V(g_tree)[v_pos_nn[i_pos]]$nclusts=tmp_data$n_points
        V(g_tree)[v_pos_nn[i_pos]]$colour=tmp_data$colour
      }
      V(g_tree)[v_pos_nn]$nclusts
      V(g_tree)[v_pos_nn]$colour
    }else{
      #take only the top 15 from data
      sub_sort=dat_clusts[1:10,]
      sub_sort$colour=top10
      
      v=match(V(g_tree)$name,as.character(sub_sort$id_Tree))
      v_pos_nn=which(!is.na(v))
      V(g_tree)[v_pos_nn]$name
      for(i_pos in 1:length(v_pos_nn))
      {
        #we need to find 
        #i_pos=2
        tmp_data=sub_sort[sub_sort$id_Tree==V(g_tree)[v_pos_nn[i_pos]]$name,]
        V(g_tree)[v_pos_nn[i_pos]]$nclusts=tmp_data$n_points
        V(g_tree)[v_pos_nn[i_pos]]$colour=tmp_data$colour
      }
      V(g_tree)[v_pos_nn]$nclusts
      V(g_tree)[v_pos_nn]$colour
      
    }
    i_level=i_level+1
    
  } #end of for(i_tmp in v_jumps)
  #Now need to add the last level of the tree
  v=match(V(g_tree)$name,"L15")
  v_pos_nn=which(!is.na(v))
  V(g_tree)[v_pos_nn]$name
  V(g_tree)[v_pos_nn]$nclusts=200000
  V(g_tree)[v_pos_nn]$colour=top10[1]
  
  v_pos_plot=which(!is.na(V(g_tree)$nclusts))
  V(g_tree)[v_pos_plot]$nclusts
  V(g_tree)[v_pos_plot]$colour
  V(g_tree)[v_pos_plot]$name
  v_x=1:(length(v_pos_plot))
  v_y=sort(V(g_tree)[v_pos_plot]$nfirms,decreasing = FALSE) 
  v_y3=v_y[1:(length(v_pos_plot)-6)]
  #v_y2=sort(V(g_plot)$nfirms,decreasing = FALSE) 
  plot(v_x,v_y)
  #now need to subset the graph and take only the ones that have values
  #--> remove nodes with NA values
  g_new=induced_subgraph(g_tree, V(g_tree)[v_pos_plot])
  #check whether is fully connected
  clusters(g_new)$no
  #take largest connected component
  gclust<-clusters(g_new, mode='weak')
  largestConnectedComponent<-induced.subgraph(g_new, V(g_new)[which(gclust$membership == which.max(gclust$csize))])
  g_plot=largestConnectedComponent
  length(V(g_plot))
  #for visualisation purposes, let us reduce the size of 6 largest ones to see the others
  V(g_plot)$size_plot=V(g_plot)$nclusts
  #for 2007 i_max=7, for 2014 i_max=5
  i_max=10 #this is for the size of the nodes
  max_nf7=sort(V(g_plot)$nclusts,decreasing = T)[i_max] 
  n_V_plot=length(V(g_plot))
  #plot(c(1:n_V_plot),V(g_plot)$nfirms)
  #plot(c(1:n_V_plot),log(V(g_plot)$nfirms))
  for(i in 1:n_V_plot)
  {
    if(V(g_plot)[i]$nclusts>max_nf7)
    {
      print(paste0(V(g_plot)[i]$name," and ",V(g_plot)[i]$nclusts))
      if(V(g_plot)[i]$name==V(g_plot)[n_V_plot]$name)
      {
        V(g_plot)[i]$size_plot=(V(g_plot)[i]$nclusts)*100000
        print(paste0(V(g_plot)[i]$name," and ",V(g_plot)[i]$nclusts))
      }else{
        V(g_plot)[i]$size_plot=(V(g_plot)[i]$nclusts)*500
      }
      print(paste0(" but now ",V(g_plot)[i]$size_plot))
    }
  }
  #sort(log(V(g_plot)$nfirms),decreasing = T)
  #y=sort((V(g_plot)$size_plot)/1200,decreasing = T)
  #plot(c(1:n_V_plot),y)
  #plot(g_plot,vertex.size=1,layout=layout_as_tree(g_plot,root = g_parent),vertex.label=NA)
  
  v_s=c(3,seq(from=2, to=0.7, by=-0.1))
  length(v_s)
  file_plot<- paste(path_plots,'/hier_tree_London.png',sep="")
  png(file_plot,height=850,width=1000)
  plot(g_plot,vertex.size=log((V(g_plot)$size_plot))/2.2,vertex.label=NA,
       layout=layout_as_tree(g_plot,root = g_parent),
       vertex.color=V(g_plot)$colour)
  
  #plot(g_plot,vertex.size=log((V(g_plot)$nfirms))/2,vertex.label=NA,layout=layout_as_tree(g_plot,root = g_parent),vertex.color=V(g_plot)$colour)
  #plot(g_plot,vertex.size=10,vertex.label=NA,layout=layout_as_tree(g_plot,root = g_parent),vertex.color=V(g_plot)$colour)
  title("Hierarchical tree, London ",cex.main=2)
  #legend(title="rank size",'topleft',legend=c(1:15),col=top15,horiz = F,pch = 19, pt.cex = v_s,bty="n",title.adj = 1.9)
  legend(-1.3,1.2,title="rank size",legend=c(1:10),col=top10,horiz = F,pch = 19, pt.cex = v_s,bty="n",
         title.adj = 1.9,cex = 1.4) #cex is used to scale size of legend
  dev.off()
  
  