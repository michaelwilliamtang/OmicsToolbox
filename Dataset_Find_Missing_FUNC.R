# finds all missing timepoints for each participant in each fiber in @dataset, output is 1 tsv files per fiber
# for mix fiber, WashoutFinal is not expected, so not marked as missing

# @dataset        e.g. genef, pcl, metaphlan, cytokine, clinical, lipids, proteomics, metabolomics
# @file_prefix    choose which version of data we use
# @wf             whether to include WashoutFinal
# @miss_sum       whether to keep track of counts of missing data points per fiber
# @source_dir     data source directory

dataset_find_missing <- function(dataset, file_prefix = "Tidy_Full",
                                 wf = T, miss_sum = T, source_dir = "Data") {
  require(tidyverse)
  
  # ranvars
  fibers = scan(file.path("Metadata", "Fibers.tsv"), character(), quote = '', sep = "\t", quiet = T)
  if (wf) {
    timepoints <- scan(file.path("Metadata", "Timepoints.tsv"), character(), quote = '', sep = "\t", quiet = T)
  } else {
    timepoints <- scan(file.path("Metadata", "Timepoints_NoWF.tsv"), character(), quote = '', sep = "\t", quiet = T)
  }
  mix_timepoints <- scan(file.path("Metadata", "Timepoints_NoWF.tsv"), character(), quote = '', sep = "\t", quiet = T)
  ids <- scan(file.path("Metadata", "Ids.tsv"), character(), quote = '', sep = "\t", quiet = T)
  
  # setting up
  if (!wf) {
    data_dir <- paste(dataset, file_prefix, "Missing_NoWF", sep = "_")
  } else {
    data_dir <- paste(dataset, file_prefix, "Missing", sep = "_")
  }
  data_dir <- file.path("Metadata", "Missing", data_dir)
  if (!dir.exists(data_dir)) {
    dir.create(data_dir)
  }
  
  # sum stuff
  if (miss_sum) {
    sum_loc <- file.path(data_dir, paste(dataset, "Missing_Sums.tsv", sep = "_"))
  }
  counts <- rep(0, length(fibers))
  names(counts) <- fibers
  
  for (fiber in fibers) {
    # exclude WAFinal as a missing in Mix
    if (fiber != "Mix") {
      tpoints <- timepoints
    } else {
      tpoints <- mix_timepoints
    }
    
    # load and setup
    file_loc <- file.path(data_dir, paste(fiber, dataset, "Missing.tsv", sep = "_"))
    load(file.path(source_dir, file_prefix, paste(file_prefix, fiber, dataset, "df.RData", sep = "_")))
    
    # setup df
    missing_df <- tibble(.rows = length(ids))
    missing_df$id <- ids
    missing_df$miss_points <- ""
    
    # subset df, find missing, write (if id not in df, mark as all missing)
    comp_df <- tidy_df %>% select(id, point) %>% unique()
    i <- 1
    for (sel_id in ids) {
      id_df <- comp_df %>% filter(id == sel_id)
      if (nrow(id_df) == 0) {
        missing_df$miss_points[i] <- paste(tpoints, collapse = "\t")
        counts[fiber] <- counts[fiber] + length(tpoints)
      } else {
        #missing_df[[sel_id]] <- paste(tpoints[which(!(tpoints %in% id_df$point))], collapse = "\t")
        missing_vec <- tpoints[which(!(tpoints %in% id_df$point))]
        missing_df$miss_points[i] <- paste(missing_vec, collapse = "\t")
      }
      i <- i + 1
      counts[fiber] <- counts[fiber] + length(missing_vec)
    }
    write_tsv(missing_df, file_loc)
  }
  
  # writing sum stuff
  if (miss_sum) {
    str_vec <- rep(0, length(fibers))
    names(str_vec) <- fibers
    for (fiber in fibers) {
      str_vec[fiber] <- paste(fiber, counts[fiber], sep = "\t")
    }
    write_lines(str_vec, sum_loc, sep = "\n")
  }
}