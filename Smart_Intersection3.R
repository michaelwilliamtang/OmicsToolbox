# a graph intersection that preserves color, checking for color collision
# this version handles collisions more naturally using edge weights, giving a positive/negative color
#   based on mean weight, but shaded to mark as a collision (red to orange, green to lightgreen)
# this version is also optimized with more vectorized color collision

smart_intersection <- function(graph1, graph2) {
  require("igraph")
  graph_new <- igraph::intersection(graph1, graph2)
  
  # optional: get rid of non-correlated nodes for easier visualization
  graph_new <- igraph::delete.vertices(graph_new, degree(graph_new) == 0)
  
  # here we avoid 0-edge graphs/0-vertex graphs that throw an error
  if (ecount(graph_new) == 0) return(graph_new)
  
  # create color columns
  E(graph_new)$color <- E(graph_new)$color_1 # assume we have edge colors, otherwise don't call this function
  if (!is.null(V(graph_new)$color_1)) { # to generalize to graphs without vertex colors
    V(graph_new)$color <- V(graph_new)$color_1 # we assume vertex colors can never be contradictory
  }
  
  # reconcile weights as mean
  if (!is.null(E(graph_new)$weight_1) && !is.null(E(graph_new)$weight_2)) {
    E(graph_new)$weight <- (E(graph_new)$weight_1 + E(graph_new)$weight_2) / 2
  }
  
  # check for color collision
  # here we assume only 2 colors, red and green, are used as base
  collisions <- which(!E(graph_new)$color_1 == E(graph_new)$color_2)
  for (e in collisions) {
    # undo abs
    weight_1 <- E(graph_new)[[e]]$weight_1
    if (E(graph_new)[[e]]$color_1 == "red") {
      weight_1 <- weight_1 * -1
    }
    weight_2 <- E(graph_new)[[e]]$weight_2
    if (E(graph_new)[[e]]$color_2 == "red") {
      weight_2 <- weight_2 * -1
    }
    
    # find weight
    weight_new <- (weight_1 + weight_2) / 2
    E(graph_new)[[e]]$weight <- abs(weight_new)
    
    if (weight_new < 0) {
      E(graph_new)[[e]]$color <- "orange"
    } else if (weight_new > 0) {
      E(graph_new)[[e]]$color <- "lightgreen"
    } else { # exactly 0 somehow
      E(graph_new)[[e]]$color <- "black"
    }
    
    # print(E(graph_new)[[e]])
  }
  print(E(graph_new)[[collisions]])
  
  return(graph_new)
}