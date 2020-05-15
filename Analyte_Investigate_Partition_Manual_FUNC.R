# investigate ONE analy as with selected partitioning type
# responder val calculated with mean of 10, 20, 30

# @dataset                e.g. genef, pcl, metaphlan, cytokine, clinical, lipids, proteomics, metabolomics
# @selected               analyte(s), default = all analytes
# @partition_dir          location of pre-computed manual partition dir
# @norm                   whether normalized with all baselines @ 0 (precomputed)

# @faceted                faceted instead of overlayed ids
# @filled                 "fill in" missing baselines with avg of present baselines (precomputed)
# @responder_label        artificially create binary division, which partition label to be treated as responder
# @partition_method       partitioning methods, current options: transparency, aggregate
# @partition_type         how to partition, current options: N_Partition, Top_N_Partition
# @only                   only include these ids
# @without                exclude these ids
# @omit_x_axis            often used for faceted graphs, omit x-axis labels for clean look
# @omit_units             often used for index or ratio analytes
# @fibers                 which fibers to graph, default = all
# @file_prefix            file prefix for data location (folder name and data file prefix)
# @graph_dir              directory for graphs
# @responder_partition    specify this partition has responders, only available if specific @partition
# @label_responders       use special labeling to single out responders, only available if @responder_partition and
#                           specified @partition
# @desc                   for documentation, should be used with any specification of partition, only, without, 
#                           multiple selected

analyte_investigate_partition_manual <- function(dataset, selected = all_analytes, partition_dir, norm = T,
                                            faceted = F, filled = F, responder_label = "responder",
                                            partition_method = "transparency", partition_type = "Manual_Partition",
                                            only = ids, without = c(), omit_x_axis = F, omit_units = F,
                                            fibers = all_fibers, file_prefix = "Tidy_Full", graph_dir = NA,
                                            responder_partition = T, label_responders = T, desc = "") {
  
  all_fibers = scan(file.path("Metadata", "Fibers.tsv"), character(), quote = '', sep = "\t", quiet = T)
  ids = scan(file.path("Metadata", "Ids.tsv"), character(), quote = '', sep = "\t", quiet = T)
  load(file.path("Data", file_prefix, paste(file_prefix, fibers[1], dataset, "df.RData", sep = "_")))
  all_analytes <- tidy_df$analyte %>% unique()
  comb_only = F # partition has no effect on combined graphs
  
  # delegate adding "Filled" to file_prefix to subfunctions
  
  for (analy in selected) {
    
    # run each fiber with its partition for that analyte; safer to specify params in case order changes
    for (fiber in fibers) {
      partition_loc <- file.path(partition_dir, paste(fiber, "_", analy, "_", dataset, ".csv", sep = ""))
      
      # consolidate in one dir
      if (is.na(graph_dir)) {
        dir_source <- "Graphs"
        dir_name <-  paste("Resp", partition_type, sep = "_")
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
      analyte_investigate_overlay(dataset = dataset, selected = analy, norm = norm, 
                          comb_only = comb_only, faceted = faceted, filled = filled, partition = partition_loc,
                          responder_label = responder_label, partition_method = partition_method,
                          only = only, without = without, omit_x_axis = omit_x_axis, omit_units = omit_units,
                          fibers = fiber, file_prefix = file_prefix, graph_dir = graph_dir,
                          responder_partition = responder_partition, label_responders = label_responders, desc = desc)
    }
  }
}