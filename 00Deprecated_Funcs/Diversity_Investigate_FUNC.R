# investigates ONE analyte in diversity data makes a fiber-comb graph and individual graphs
#   per fiber showing the participants
# updated with units, takes care of protein full names and clinical units in particular

# @analy    which analyte (shannon_diversity or simpson_diversity) you want
# @norm     whether normalized with all baselines @ 0
# @fibers                 which fibers to graph, default = all

diversity_investigate <- function(analy, norm, fibers) {

  require(tidyverse)
  require(plotrix)
  
  fibers = c("Arabinoxylan","LCInulin","Mix")
  file_prefix <- "Tidy_Species_Diversity"
  data_dir <- file.path("Data", file_prefix)
  dataset <- "metaphlan"
  
  dir_source <- "Graphs"
  dir_name <- "Diversity_Analyses"
  if (norm) dir_name <- paste(dir_name, "Norm", sep = "_")
  # if (faceted) dir_name <- paste(dir_name, "Faceted", sep = "_")
  if (length(fibers) == 1) dir_name <- paste(dir_name, fibers, sep = "_")
  
  graph_dir <- paste(dir_name, analy, sep = "_")
  if (!dir.exists(graph_dir))
    dir.create(graph_dir)
  
  comb_df <- data.frame()
  
  analy_full <- analy
  unit <- ""
  
  print(paste(analy_full, unit, sep = ""))
  
  for (fiber in fibers) {
    print(fiber)
    
    # loading data
    load(file.path(data_dir, paste("Tidy", "species_diversity", fiber, dataset, "df.RData", sep = "_")))
    tidy_cyto <- tidy_df %>% filter(analyte == analy)
    
    # graphing ids together, averaged
    if (norm) {
      tmp_cyto <- tidy_cyto %>%
        group_by(point) %>%
        dplyr::summarise(mean_parts = mean(renorm_val), std_error = std.error(val)) %>% # std error cannot be based on norm
        ungroup()
    } else {
      tmp_cyto <- tidy_cyto %>%
        group_by(point) %>%
        dplyr::summarise(mean_parts = mean(val), std_error = std.error(val)) %>% # std error cannot be based on norm
        ungroup()
    }
    
    # saving dfs
    tmp_cyto <- cbind(tmp_cyto, "fiber" = rep(fiber, nrow(tmp_cyto)))
    comb_df <- rbind(comb_df, tmp_cyto)
    
    # graphing ids separately
    if (norm) {
      pdf(file.path(graph_dir, paste(analy, dataset, fiber, "Ids.pdf", sep = "_")), width = 9, height = 6)
      print(tidy_cyto %>%
        ggplot(aes(x = point, y = renorm_val, group = id, color = id)) +
        # geom_point() +
        geom_line() +
        theme(
          axis.text.x = element_text(size=10, angle=45, hjust = 1)) +
        ylab(paste(analy_full, unit, sep = "")))
      dev.off()
    } else {
      pdf(file.path(graph_dir, paste(analy, dataset, fiber, "Ids.pdf", sep = "_")), width = 9, height = 6)
      print(tidy_cyto %>%
              ggplot(aes(x = point, y = val, group = id, color = id)) +
              # geom_point() +
              geom_line() +
              theme(
                axis.text.x = element_text(size=10, angle=45, hjust = 1)) +
              ylab(paste(analy_full, unit, sep = "")))
      dev.off()
    }
  }
  
  if (length(fibers) == 1) return("Only 1 fiber, so no comb graph")
  
  # with error bars, with comb_df and working dodging
  pd <- position_dodge(width = 0.2)
  
  pdf(file.path(graph_dir, paste(analy, "_", dataset, "Comb_Avg.pdf", sep = "")), width = 9, height = 6)
  
  # normalizing with each fiber baseline
  # assumed sorted, so baseline is first timepoint of each analyte-fiber set, and fibers adjacent
  
  print(comb_df %>%
    ggplot(position = pd, aes(x = point, y = mean_parts, group = fiber, color = fiber)) +
    geom_line(position = pd) +
    geom_errorbar(aes(x = point, ymin = mean_parts - std_error, ymax = mean_parts + std_error),
                  width = .5, position = pd) +
    theme(
      axis.text.x = element_text(size=10, angle=45, hjust = 1)) +
    ylab(paste(analy_full, unit, sep = "")))
  dev.off()
}