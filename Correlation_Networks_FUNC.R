# correlation networks for 2 omics/datasets, as well as graph algebra permutations, from given correlations
# graph types: main (shows all correlations), only (only shows cross-omic edges), all (shows
#   all edges between the vertices of cross-omic correlations); up to 30 graphs total
# non-correlating vertices omitted in all cases

# @dataset1/2             e.g. genef, pcl, metaphlan, cytokine, clinical, lipids, proteomics, 
#                         metabolomics
# @fibers                 fibers to include, default = all
# @source_dir             data location
# @file_prefix            prefix of file name
# @file_postfix           postfix of file name, optional
# @desc                   optional, for documentation, needs to describe clearly

correlation_networks <- function(dataset1, dataset2, fibers = all_fibers, 
                             source_dir = "Data", file_prefix = "Tidy_Full", file_postfix = "", desc = "") {
  
  require(Mfuzz)
  require(matrixStats)
  require(icesTAF)
  require(Hmisc)
  require(igraph)
  require(numform)
  require(tidyverse)
  # require(impute)
  
  all_fibers <- scan(file.path("Metadata", "Fibers.tsv"), character(), quote = '', sep = "\t", quiet = T)
  file_prefix1 <- "Tidy_Full" # aux directory to get analyte lists
  load(file.path("Data", file_prefix1, paste(file_prefix1, fibers[1], dataset1, "df.RData", sep = "_")))
  all_analytes1 <- tidy_df$analyte %>% unique()
  load(file.path("Data", file_prefix1, paste(file_prefix1, fibers[1], dataset2, "df.RData", sep = "_")))
  all_analytes2 <- tidy_df$analyte %>% unique()
  
  # creating dirs
  fiber_dir <- file.path("Graphs", paste("Correlation_Networks", dataset1, dataset2,
                                        desc, sep = "_")) # quick dir name change
  fiber_dir1 <- file.path(fiber_dir, "Correlations") # quick dir name change
  fiber_dir2 <- file.path(fiber_dir, "Comparisons") # quick dir name change
  if (!dir.exists(fiber_dir))
    dir.create(fiber_dir)
  if (!dir.exists(fiber_dir1))
    dir.create(fiber_dir1)
  if (!dir.exists(fiber_dir2))
    dir.create(fiber_dir2)
  
  graph_list <- list() # for comparison graph storage
  
  for (fiber in fibers) {
    print(fiber) # verbose
    
    # scan correlation data
    file_name <- paste(file_prefix, fiber, sep = "_")
    if (length(file_postfix) > 0) file_name <- paste(file_name, file_postfix, sep = "_")
    cor.data <- read.table(file.path(source_dir, file_name), sep = "\t")
    cor.data <- data.matrix(cor.data)
    
    ##################################################################
    ### making network
    ##################################################################
    
    network <- graph.adjacency(abs(cor.data), weighted=T, mode="undirected", diag=F)
    if (ecount(network) == 0) next
    
    # corr to edge colors
    edges <- data.frame(get.edgelist(network))
    edges <- transform(edges, X1 = as.character(X1), X2 = as.character(X2)) # currently factor, which we don't want
    correlation_values <- rep(0, times = nrow(edges)) # allocate instead of growing for efficiency
    for (i in 1:nrow(edges)){
      correlation_values[[i]] <- cor.data[edges[i,2],edges[i,1]]
    }
    E(network)$color <- ifelse(correlation_values < 0, "red", "blue")
    # E(network)$weight <- correlation_values
    
    # vertex colors
    num_dataset1 <- 130 #length(which(rownames(cor.data) %in% all_analytes1))
    num_dataset2 <- 107 #length(which(rownames(cor.data) %in% all_analytes2))
    if (num_dataset1 == 0) next # catching no-analy1 situation
    if (num_dataset2 == 0) next # catching no-analy2 situ
    colors <- c(rep("#0006f5ff", length = num_dataset1), # "#0006f5ff"
                rep("yellow", length = num_dataset2))
    # colors <- c(rep("lightblue", length = vcount(network))) # only 1 dataset
    V(network)$color <- colors
    
    # # removing vertices that are not connected to a clinical
    # clin_verts <- V(network)[1:length(clinical_obs)]
    # non_clin_edges <- E(network)[!.inc(clin_verts)]
    # network2 <- delete_edges(network, non_clin_edges)
    # network2 <- delete_vertices(network2, degree(network2) == 0)
    # 
    # network <- delete_vertices(network, degree(network) == 0)
    
    network <- igraph::delete.vertices(network, degree(network) == 0)
    # this used to be length(clinical_obs) but we deleted unconnected verts
    # clin_num <- max(which(V(network)$color == "yellow"))
    
    ##################################################################
    ### plotting full
    ##################################################################
    
    pdf(file.path(fiber_dir1, paste(fiber, "Correlations.pdf", sep = "_")))
    
    ly <- layout.fruchterman.reingold(network, dim=2, grid="nogrid")
    par(bg="white", mar=c(0,0,0,0))
    set.seed(4)
    plot(network,
         vertex.size=7,
         # vertex.color=V(network)$color,
         # vertex.label = NA,
         vertex.label.cex=0.5,
         vertex.label.color="black",
         vertex.frame.color="black",
         layout = ly,
         edge.width = E(network)$weight * 5
    )
    dev.off()
    
    # storage
    graph_list[[fiber]] <- network
    
    ##################################################################
    ### getting and plotting interesting edges (only)
    ##################################################################
    
    # # getting list of clinical-protein edges
    # clin_verts <- V(network)[1 : clin_num]
    # clin_edges <- E(network)[.inc(clin_verts)]
    # protein_verts <- V(network)[(clin_num + 1) : vcount(network)]
    # both_edges <- clin_edges[.inc(protein_verts)]
    # 
    # if (length(both_edges) == 0) next # catching no-both situation
    # 
    # rm_edges <- E(network)[-both_edges]
    # network_tmp <- delete_edges(network, rm_edges)
    # 
    # # tidying
    # network_edges <- cbind(as_edgelist(network_tmp), edge_attr(network_tmp, "color"))
    # colnames(network_edges) <- c("node_1", "node_2", "color")
    # 
    # # writing
    # f <- file.path(fiber_dir1, paste(fiber, "Partial_Edges.csv", sep = "_"))
    # write.csv(network_edges, file = f, row.names = F)
    # 
    # # skeleton plot
    # pdf(file.path(fiber_dir1, paste(fiber, "Partial_Correlations_Only.pdf", sep = "_")))
    # network_tmp <- igraph::delete.vertices(network_tmp, degree(network_tmp) == 0)
    # ly_tmp <- layout.fruchterman.reingold(network_tmp, dim=2, grid="nogrid")
    # par(bg="white", mar=c(0,0,0,0))
    # set.seed(4)
    # plot(network_tmp,
    #      vertex.size=7,
    #      # vertex.color=V(network)$color,
    #      vertex.label.cex=0.5,
    #      vertex.label.color="black",
    #      vertex.frame.color="black",
    #      layout = ly_tmp,
    #      edge.width = E(network_tmp)$weight * 5
    # )
    # dev.off()
    
    ##################################################################
    ### plotting partial all edges
    ##################################################################
    
    # both_verts <- ends(network, both_edges, names = F) # nums faster than names
    # both_verts <- both_verts %>% as.vector() %>% unique()
    # 
    # network_tmp3 <- induced_subgraph(network, both_verts)
    # 
    # # plotting
    # pdf(file.path(fiber_dir1, paste(fiber, "Partial_Correlations_All.pdf", sep = "_")))
    # 
    # # network_tmp3 <- igraph::delete.vertices(network_tmp3, degree(network_tmp3) == 0)
    # ly_tmp3 <- layout.fruchterman.reingold(network_tmp3, dim=2, grid="nogrid")
    # par(bg="white", mar=c(0,0,0,0))
    # set.seed(4)
    # plot(network_tmp3,
    #      vertex.size=7,
    #      # vertex.color=V(network)$color,
    #      vertex.label.cex=0.5,
    #      vertex.label.color="black",
    #      vertex.frame.color="black",
    #      layout = ly_tmp3,
    #      edge.width = E(network_tmp3)$weight * 5
    # )
    # dev.off()
    
    ##################################################################
    ### plotting partial cluster
    ##################################################################
    
    # # getting clusters with clins
    # clu <- clusters(network)
    # 
    # # clin_clu <- rep(NA, clu$no) # allocating
    # # clin_clu <- clin_clu[!is.na(clin_clu)] # getting rid of allocated unused spaces
    # 
    # # getting selected nodes
    # clin_clu <- unique(clu$membership[1:clin_num])
    # prot_clu <- unique(clu$membership[(clin_num + 1):vcount(network)])
    # both_clu <- base::intersect(clin_clu, prot_clu)
    # # sel_clu <- names(clu$membership)[which(clu$membership %in% clin_clu)] # not preferred; vert names
    # sel_clu <- which(clu$membership %in% both_clu) # actually prefer nums
    # 
    # network_tmp2 <- induced_subgraph(network, sel_clu)
    # 
    # # plotting
    # pdf(file.path(fiber_dir1, paste(fiber, "Partial_Correlations_Cluster.pdf", sep = "_")))
    # 
    # # network_tmp2 <- igraph::delete.vertices(network_tmp2, degree(network_tmp2) == 0)
    # ly_tmp2 <- layout.fruchterman.reingold(network_tmp2, dim=2, grid="nogrid")
    # par(bg="white", mar=c(0,0,0,0))
    # set.seed(4)
    # plot(network_tmp2,
    #      vertex.size=7,
    #      # vertex.color=V(network)$color,
    #      vertex.label.cex=0.5,
    #      vertex.label.color="black",
    #      vertex.frame.color="black",
    #      layout = ly_tmp2
    # )
    # dev.off()
  }
  
  ##################################################################
  ### comparisons
  ##################################################################
  
  source("Smart_Intersection3.R")
  # these are capital snake case for better pdf printing later
  comparisons <- list (
    Intersect_All = smart_intersection(smart_intersection(graph_list[["LCInulin"]], graph_list[["Mix"]]), graph_list[["Arabinoxylan"]]),
    Unique_Ara = difference(difference(graph_list[["Arabinoxylan"]], graph_list[["LCInulin"]]), graph_list[["Mix"]]),
    Unique_Inu = difference(difference(graph_list[["LCInulin"]], graph_list[["Arabinoxylan"]]), graph_list[["Mix"]]),
    Unique_Mix = difference(difference(graph_list[["Mix"]], graph_list[["LCInulin"]]), graph_list[["Arabinoxylan"]]),
    Unique_Arainu = difference(smart_intersection(graph_list[["LCInulin"]], graph_list[["Arabinoxylan"]]), graph_list[["Mix"]]),
    Unique_Aramix = difference(smart_intersection(graph_list[["Mix"]], graph_list[["Arabinoxylan"]]), graph_list[["LCInulin"]]),
    Unique_Inumix = difference(smart_intersection(graph_list[["LCInulin"]], graph_list[["Mix"]]), graph_list[["Arabinoxylan"]])
  )
  
  # check to make sure comparisons sum to total unique edges (optional)
  count <- 0
  for (c in comparisons) {
    count <- count + length(E(c))
  }
  test_union <- graph.union(graph.union(graph_list[["LCInulin"]], graph_list[["Mix"]]), graph_list[["Arabinoxylan"]])
  count2 <- length(E(test_union))
  print(count == count2)
  
  
  # plotting
  for (c in names(comparisons)) {
    
    network <- comparisons[[c]] # convenience and compatibility
    
    # here we avoid 0-edge graphs which will become 0-vertex graphs and throw an error
    if (ecount(network) == 0) next
    
    par(bg="white", mar=c(0,0,0,0))
    
    # optional: get rid of non-correlated nodes for easier visualization; right now it is also done in smart_intersection
    network <- igraph::delete.vertices(network, degree(network) == 0)
    # network <- delete.vertices(simplify(network), degree(simplify(network)) == 0)
    
    # plot
    pdf(file.path(fiber_dir2, paste(c, "Comparison.pdf", sep = "_")))
    ly <- layout.fruchterman.reingold(network,dim=2,grid="nogrid")
    set.seed(4)
    plot(network,
         vertex.size=7,
         vertex.label.cex=0.5,
         vertex.label.color="black",
         vertex.frame.color="black",
         #vertex.label = NA,
         layout = ly,
         edge.width = E(network)$weight * 5
    )
    
    dev.off()
    
    # no guarentee the comp graph exists
    # if (length(which(V(network)$color == "yellow")) == 0) next # catching no-clinical situation
    # clin_num <- max(which(V(network)$color == "yellow"))
    
    ##################################################################
    ### getting interesting edges
    ##################################################################
    
    # clin_verts <- V(network)[1 : clin_num]
    # clin_edges <- E(network)[.inc(clin_verts)]
    # protein_verts <- V(network)[(clin_num + 1) : vcount(network)]
    # both_edges <- clin_edges[.inc(protein_verts)]
    # 
    # if (length(both_edges) == 0) next # catching no-both situation
    # 
    # rm_edges <- E(network)[-both_edges]
    # network_tmp <- delete_edges(network, rm_edges)
    # 
    # # tidying
    # network_edges <- cbind(as_edgelist(network_tmp), edge_attr(network_tmp, "color"))
    # colnames(network_edges) <- c("node_1", "node_2", "color")
    # 
    # # writing
    # f <- file.path(fiber_dir2, paste(c, "Partial_Edges.csv", sep = "_"))
    # write.csv(network_edges, file = f, row.names = F)
    
    ##################################################################
    ### plotting interesting edges
    ##################################################################
    
    # pdf(file.path(fiber_dir2, paste(c, "Partial_Comparison_Only.pdf", sep = "_")))
    # 
    # network_tmp <- igraph::delete.vertices(network_tmp, degree(network_tmp) == 0)
    # ly_tmp <- layout.fruchterman.reingold(network_tmp, dim=2, grid="nogrid")
    # par(bg="white", mar=c(0,0,0,0))
    # set.seed(4)
    # plot(network_tmp,
    #      vertex.size=7,
    #      # vertex.color=V(network)$color,
    #      vertex.label.cex=0.5,
    #      vertex.label.color="black",
    #      vertex.frame.color="black",
    #      layout = ly_tmp,
    #      edge.width = E(network_tmp)$weight * 5
    # )
    # 
    # dev.off()
    # 
    ##################################################################
    ### plotting partial all edges
    ##################################################################
    
    # both_verts <- ends(network, both_edges, names = F) # nums faster than names
    # both_verts <- both_verts %>% as.vector() %>% unique()
    # 
    # network_tmp3 <- induced_subgraph(network, both_verts)
    # 
    # # plotting
    # pdf(file.path(fiber_dir2, paste(c, "Partial_Comparison_All.pdf", sep = "_")))
    # 
    # # network_tmp3 <- igraph::delete.vertices(network_tmp3, degree(network_tmp3) == 0)
    # ly_tmp3 <- layout.fruchterman.reingold(network_tmp3, dim=2, grid="nogrid")
    # par(bg="white", mar=c(0,0,0,0))
    # set.seed(4)
    # plot(network_tmp3,
    #      vertex.size=7,
    #      # vertex.color=V(network)$color,
    #      vertex.label.cex=0.5,
    #      vertex.label.color="black",
    #      vertex.frame.color="black",
    #      layout = ly_tmp3,
    #      edge.width = E(network_tmp3)$weight * 5
    # )
    # dev.off()
    
    ##################################################################
    ### plotting partial cluster
    ##################################################################
    
    # # getting clusters with clins
    # clu <- clusters(network)
    # 
    # # clin_clu <- rep(NA, clu$no) # allocating
    # # clin_clu <- clin_clu[!is.na(clin_clu)] # getting rid of allocated unused spaces
    # 
    # # getting selected nodes
    # clin_clu <- unique(clu$membership[1:clin_num])
    # prot_clu <- unique(clu$membership[(clin_num + 1):vcount(network)])
    # both_clu <- base::intersect(clin_clu, prot_clu)
    # # sel_clu <- names(clu$membership)[which(clu$membership %in% clin_clu)] # not preferred; vert names
    # sel_clu <- which(clu$membership %in% both_clu) # actually prefer nums
    # 
    # network_tmp2 <- induced_subgraph(network, sel_clu)
    # 
    # # plotting
    # pdf(file.path(fiber_dir2, paste(c, "Partial_Comparison_Cluster.pdf", sep = "_")))
    # 
    # # network_tmp2 <- igraph::delete.vertices(network_tmp2, degree(network_tmp2) == 0)
    # ly_tmp2 <- layout.fruchterman.reingold(network_tmp2, dim=2, grid="nogrid")
    # par(bg="white", mar=c(0,0,0,0))
    # set.seed(4)
    # plot(network_tmp2,
    #      vertex.size=7,
    #      # vertex.color=V(network)$color,
    #      vertex.label.cex=0.5,
    #      vertex.label.color="black",
    #      vertex.frame.color="black",
    #      layout = ly_tmp2
    # )
    # dev.off()
    
  }
}