# investigate ONE analy as 3 groups on 1 graph, with top, middle, bottom groups of responders
# responder val calculated with mean of 10, 20, 30

# @dataset      e.g. genef, pcl, metaphlan, cytokine, clinical, lipids, proteomics, metabolomics
# @analy        desired analyte in the dataset
# @norm         whether normalized with all baselines @ 0 (precomputed)
# @overwrite    whether to redo partition if already exists
# @filled       "fill in" missing baselines with avg of present baselines (precomputed)
# @only         only include these ids
# @without      exclude these ids

responder_tripartition <- function(dataset, analy, norm = T, overwrite = F, fibers = all_fibers,
                                   filled = F, only = ids, without = c()) {
  require(tidyverse)
  require(plotrix)
  
  if (filled) file_prefix <- "Tidy_Full_Filled"
  else  file_prefix <- "Tidy_Full"
  
  all_fibers = scan(file.path("Metadata", "Fibers.tsv"), character(), quote = '', sep = "\t", quiet = T)
  ids = scan(file.path("Metadata", "Ids.tsv"), character(), quote = '', sep = "\t", quiet = T)
  
  partition_name <- paste("Partition", analy, dataset, sep = "_")
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
    
    # making groups
    response <- factor(levels = c("top", "middle", "bottom"))
    resp_df$partition <- levels(response)[(0:(nrow(resp_df) - 1) %/% (nrow(resp_df)/ 3)) + 1] # splitting into 3 groups
    
    # save
    write.csv(resp_df, partition_loc)
  }
  return(partition_dir)
}