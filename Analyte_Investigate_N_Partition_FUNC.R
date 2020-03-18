# investigate ONE analy as 3 groups on 1 graph, with top, middle, bottom groups of responders
# responder val calculated with mean of 10, 20, 30

# @dataset                e.g. genef, pcl, metaphlan, cytokine, clinical, lipids, proteomics, metabolomics
# @selected               analyte(s), default = all analytes
# @norm                   whether normalized with all baselines @ 0 (precomputed)
# @N                      number of partitions
# @overwrite              whether to redo partition if already exists

# @faceted                faceted instead of overlayed ids
# @filled                 "fill in" missing baselines with avg of present baselines (precomputed)
# @responder_label        artificially create binary division, which partition label to be treated as responder
# @only                   only include these ids
# @without                exclude these ids
# @omit_x_axis            often used for faceted graphs, omit x-axis labels for clean look
# @fibers                 which fibers to graph, default = all
# @graph_dir              directory for graphs
# @responder_partition    specify this partition has responders, only available if specific @partition
# @label_responders       use special labeling to single out responders, only available if @responder_partition and
#                           specified @partition
# @desc                   for documentation, should be used with any specification of partition, only, without, 
#                           multiple selected

analyte_investigate_N_partition <- function(dataset, selected = all_analytes, norm = T, N = 3, overwrite = F,
                                            faceted = F, filled = F, responder_label = "top",
                                            only = ids, without = c(), omit_x_axis = F,
                                            fibers = all_fibers, graph_dir = NA,
                                            responder_partition = T, label_responders = T, desc = "") {
  
  all_fibers = scan(file.path("Metadata", "Fibers.tsv"), character(), quote = '', sep = "\t", quiet = T)
  ids = scan(file.path("Metadata", "Ids.tsv"), character(), quote = '', sep = "\t", quiet = T)
  load(file.path("Data", prefix, paste(prefix, fibers[1], dataset, "df.RData", sep = "_")))
  all_analytes <- tidy_df$analyte %>% unique()
  comb_only = F # partition has no effect on combined graphs
  
  for (analy in selected) {
    # get partitions; safer to specify params in case order changes
    partition_dir <- N_partition(dataset = dataset, analy = analy, norm = norm, N = N, overwrite = overwrite,
                                            fibers = fibers, filled = filled, only = only, without = without)
    
    
    # run each fiber with its partition for that analyte; safer to specify params in case order changes
    for (fiber in fibers) {
      partition_loc <- file.path(partition_dir, paste(fiber, "_", analy, "_", dataset, ".csv", sep = ""))
      
      # consolidate in one dir
      if (is.na(graph_dir)) {
        dir_source <- "Graphs"
        dir_name <-  paste("Resp", N, "Partition", sep = "_")
        if (norm) dir_name <- paste(dir_name, "Norm", sep = "_")
        if (faceted) dir_name <- paste(dir_name, "Faceted", sep = "_")
        if (filled) dir_name <- paste(dir_name, "Filled", sep = "_")
        if (length(selected) == 1) dir_name <- paste(dir_name, selected, sep = "_")
        if (length(fibers) == 1) dir_name <- paste(dir_name, fibers, sep = "_")
        dir_name <- paste(dir_name, dataset, sep = "_")
        if (desc != "") dir_name <- paste(dir_name, desc, sep = "_")
        graph_dir <- file.path(dir_source, dir_name)
      }
      if (!dir.exists(graph_dir))
        dir.create(graph_dir)
      
      # make graphs
      analyte_investigate(dataset = dataset, selected = analy, norm = norm, 
                          comb_only = comb_only, faceted = faceted, filled = filled, partition = partition_loc,
                          responder_label = responder_label,
                          only = only, without = without, omit_x_axis = omit_x_axis, 
                          fibers = fiber, graph_dir = graph_dir,
                          responder_partition = responder_partition, label_responders = label_responders, desc = desc)
    }
  }
}