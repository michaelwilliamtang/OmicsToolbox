# investigate ONE analy as 2 groups on 1 graph, marks top N as responder, others as nonresponder
# responder val calculated with mean of 10, 20, 30

# @dataset          e.g. genef, pcl, metaphlan, cytokine, clinical, lipids, proteomics, metabolomics
# @analy            desired analyte in the dataset
# @norm             whether normalized with all baselines @ 0 (precomputed)
# @N                number of partitions
# @overwrite        whether to redo partition if already exists
# @filled           "fill in" missing baselines with avg of present baselines (precomputed)
# @only             only include these ids
# @without          exclude these ids
# @fibers           which fibers to graph, default = all
# @file_prefix      file prefix for data location (folder name and data file prefix)

top_N_partition <- function(dataset, analy, norm = T, N = 3, overwrite = F, filled = F, only = ids, without = c(),
                        fibers = all_fibers, file_prefix = "Tidy_Full") {
  require(tidyverse)
  require(plotrix)
  
  if (filled) file_prefix <- paste(file_prefix, "Filled", sep = "_")
  
  all_fibers = scan(file.path("Metadata", "Fibers.tsv"), character(), quote = '', sep = "\t", quiet = T)
  ids = scan(file.path("Metadata", "Ids.tsv"), character(), quote = '', sep = "\t", quiet = T)
  
  partition_name <- paste("Top", N, "Partition", analy, dataset, sep = "_")
  partition_dir <- file.path("Metadata", "Partitions", partition_name)
  if (dir.exists(partition_dir)) {
    print("Partition already created")
    if (overwrite) print("Overwriting")
    else return(partition_dir)
  }
  else {
    dir.create(partition_dir)
    print("Creating new partition")
  }
  
  for (fiber in fibers) {
    
    partition_loc <- file.path(partition_dir, paste(fiber, "_", analy, "_", dataset, ".csv", sep = ""))
    
    # load
    load(file.path("Data", file_prefix, paste(file_prefix, fiber, dataset, "df.RData", sep = "_")))
    tidy_df <- tidy_df %>% filter(analyte == analy) %>% filter(id %in% only) %>% filter(!(id %in% without))
    
    # compute resp values; method = average 10, 20, 30 timepoints
    if (norm) {
      resp_df <- tidy_df %>% filter(point == "30" | point == "20" | point == "10") %>%
        select(-c(val, analyte)) %>% # necessary for clean/proper spread
        spread(point, renorm_val, sep = "_")
    } else {
      resp_df <- tidy_df %>% filter(point == "30" | point == "20" | point == "10") %>%
        select(-c(renorm_val, analyte)) %>%
        spread(point, val, sep = "_") 
    }
    resp_df <- resp_df %>%
      mutate(avg = rowMeans(dplyr::select(., point_30, point_20, point_10), na.rm = T)) %>%
      select(id, avg) %>%
      arrange(desc(avg))
    
    if (N < 1) return ("Invalid N")
    resp_df$partition <- c(rep("responder", times = N), rep("nonresponder", times = nrow(resp_df) - N))
    
    # save
    write.csv(resp_df, partition_loc)
  }
  return(partition_dir)
}