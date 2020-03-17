# investigates ONE analyte in an omic/dataset and makes a fiber-comb graph and individual graphs
#   per fiber showing the participants
# renorms centered at part means
# adds normal_range ranges for clinicals
# units, full names and clinical units in particular

# @dataset    string e.g. genef, pcl, metaphlan, cytokine, clinical, lipids, proteomics, 
#             metabolomics
# @selected   vector of analyte string(s), default = all analytes
# @renorm     whether renormalized with all baselines @ 0
# @partition  how to divide the data, default = ids separate
# @only       only include these ids, default = include all
# @without    exclude these ids, default = exclude none (can only use only OR without, not both together)
# @list_name  optional, for documentation, needs to describe clearly

# "clean-up" changes:
# removed x-axis tick labels, x-axis title (point), y-axis "baselines normalized to mean," legend, different id colors 
# (removal of "baselines norm" ends up affecting both faceted and Comb_Avg graphs)
# filled: due to large number of missing baselines, fill in baselines with avg baseline of all participants as compromise

analyte_investigate_faceted <- function(dataset, selected = all_analytes, 
                                        faceted = T, renorm = T, filled = T, partition = NA,
                                        only = ids, without = c(), list_name = "") {
  require(tidyverse)
  require(plotrix)
  
  prefix <- "Tidy_Full"
  if (filled) prefix <- paste(prefix, "Filled", sep = "_")

  fibers = scan(file.path("Metadata", "Fibers.tsv"), character(), quote = '', sep = "\t")
  ids = scan(file.path("Metadata", "Ids.tsv"), character(), quote = '', sep = "\t")
  load(file.path("Data", prefix, paste(prefix, fibers[1], dataset, "df.RData", sep = "_")))
  all_analytes <- tidy_df$analyte %>% unique()
  
  dir_source <- "Graphs"
  dir_name <- "Analyte"
  if (faceted) dir_name <- paste(dir_name, "Faceted", sep = "_")
  if (filled) dir_name <- paste(dir_name, "Filled", sep = "_")
  if (renorm) dir_name <- paste(dir_name, "Renorm", sep = "_")
  if (length(selected) == 1) {
    graph_dir <- file.path(dir_source, paste(dir_name, selected, dataset, sep = "_"))
  } else {
    graph_dir <- file.path(dir_source, paste(dir_name, dataset, list_name, sep = "_"))
  }
  if (!dir.exists(graph_dir))
    dir.create(graph_dir)
  
  if (dataset == "clinical") {
    # add units
    clin_names <- scan(file.path("Metadata", "Clin_Names_Safe.tsv"), character(), quote = '', sep = "\t")
    clin_names_full <- scan(file.path("Metadata", "Clin_Names.tsv"), character(), quote = '', sep = "\t")
    clin_units <- scan(file.path("Metadata", "Clin_Units.tsv"), character(), quote = '', sep = "\t")
    names(clin_units) <- clin_names
    names(clin_names_full) <- clin_names
    
    # also reading for normal_range ranges
    clin_mins <- scan(file.path("Metadata", "Clin_Healthy_Ranges_min.tsv"), double(), quote = '', sep = "\t")
    clin_maxs <- scan(file.path("Metadata", "Clin_Healthy_Ranges_max.tsv"), double(), quote = '', sep = "\t")
    names(clin_mins) <- clin_names
    names(clin_maxs) <- clin_names
  }
  
  if (dataset == "proteomics") {
    # loading data
    load("Data/Full/Proteomics.RData")
    proteo <- Proteomics_data_1pept_Log2_filt80perc_ImptND_Combat
    proteo$Protein.name <- as.character(proteo$Protein.name)
  }
  
  # load partition data
  if (!is.na(partition)) {
    partition_df <- read.csv(partition, stringsAsFactors = F)
    partition_vec <- partition_df$partition
    names(partition_vec) <- partition_df$participant
  }
  
  for (analy in selected) {
    
    comb_df <- data.frame()
  
    # setting units
    if (dataset == "clinical") {
      if (is.na(clin_units[analy])) {
        unit <- ""
      } else {
        unit <- paste(" (", clin_units[analy], ")", sep = "")
      }
      # if (renorm) {
      #   unit <- paste(unit, ",\nbaselines normalized to mean", sep = "")
      # }
    } else {
      if (renorm) {
        unit <- " (log2 normalized expression)"
      } else {
        unit <- " (log2 expression)"
      }
    }
    
    if (dataset == "clinical") {
      analy_full <- clin_names_full[analy]
    } else {
      analy_full <- analy
    }
    
    if (dataset == "proteomics") {
      analy_full <- proteo[gsub("[.].$", "", analy), "Protein.name"]
    }

    print(paste(analy_full, unit, sep = ""))
    
    for (fiber in fibers) {
      print(fiber)
      
      # loading data
      load(file.path("Data", prefix, paste(prefix, fiber, dataset, "df.RData", sep = "_")))
      tidy_cyto <- tidy_df %>% filter(analyte == analy) %>% filter(id %in% only) %>% filter(!(id %in% without))
      
      # graphing ids together, averaged
      if (renorm) {
        tidy_cyto <- na.omit(tidy_cyto) # for missing baselines creating NAs in renormalization
        tmp_cyto <- tidy_cyto %>%
          group_by(point) %>%
          dplyr::summarise(mean_parts = mean(renorm_val), std_error = std.error(renorm_val)) %>% # std error cannot be based on renorm
          ungroup()
      } else {
        tmp_cyto <- tidy_cyto %>%
          group_by(point) %>%
          dplyr::summarise(mean_parts = mean(val), std_error = std.error(val)) %>% # std error cannot be based on renorm
          ungroup()
      }
      
      # saving dfs
      tmp_cyto <- cbind(tmp_cyto, "fiber" = rep(fiber, nrow(tmp_cyto)))
      comb_df <- rbind(comb_df, tmp_cyto)
      
      if (dataset == "clinical") {
        # setting up background rects
        if (renorm) {
          data_min <- min(tidy_cyto$renorm_val)
          data_max <- max(tidy_cyto$renorm_val)
        } else {
          data_min <- min(tidy_cyto$val)
          data_max <- max(tidy_cyto$val)
        }
        # data_min <- data_mins[[analy]]
        # data_max <- data_maxs[[analy]]
        if (data_min >= clin_maxs[analy] || data_max <= clin_mins[analy]) { # case of no normal_range
          ranges <- data.frame(ystart = data_min,
                               yend = data_max,
                               col = "outside_normal")
        } else if (data_min >= clin_mins[analy] && data_max <= clin_maxs[analy]) { # case of normal_range only
          ranges <- data.frame(ystart = data_min,
                               yend = data_max,
                               col = "normal_range")
        } else if (data_min >= clin_mins[analy]) { # case of normal_range + outside_normal upper
          ranges <- data.frame(ystart = c(data_min, clin_maxs[analy]),
                               yend = c(clin_maxs[analy], data_max),
                               col = c("normal_range", "outside_normal"))
        } else if (data_max <= clin_maxs[analy]) { # case of normal_range + outside_normal lower
          ranges <- data.frame(ystart = c(data_min, clin_mins[analy]),
                               yend = c(clin_mins[analy], data_max),
                               col = c("outside_normal", "normal_range"))
        } else { # case normal_range + outside_normal on both sides
          ranges <- data.frame(ystart = c(data_min, clin_mins[analy], clin_maxs[analy]),
                               yend = c(clin_mins[analy], clin_maxs[analy], data_max),
                               col = c("outside_normal", "normal_range", "outside_normal"))
        }
      }
      
      # partitioning
      if (is.na(partition)) {
        tidy_cyto$partition = tidy_cyto$id
      } else {
        tidy_cyto$partition = partition_vec[tidy_cyto$id]
      }
      
      # graphing ids separately
      pdf(file.path(graph_dir, paste(analy, dataset, fiber, "Ids.pdf", sep = "_")), width = 9, height = 6)
      
      plot <- tidy_cyto %>% ggplot()
      if (dataset == "clinical") plot <- plot + 
        geom_rect(data = ranges, aes(ymin = ystart, ymax = yend, xmin = -Inf, xmax = Inf, fill = col), alpha = 0.4) +
        scale_fill_manual(values = c(outside_normal = "#ffcccb", normal_range = "#c6ff95"))
      if (renorm && faceted) plot <- plot + 
        geom_line(aes(x = point, y = renorm_val, group = id))
      else if (renorm && !faceted)  plot <- plot + 
        geom_line(aes(x = point, y = renorm_val, group = id, color = id))
      else if (!renorm && faceted) plot <- plot + 
        geom_line(aes(x = point, y = val, group = id))
      else plot <- plot + 
        geom_line(aes(x = point, y = val, group = id, color = id))
      # clean theme
      plot <- plot + 
        theme(
          # axis.text.x = element_text(size=10, angle=45, hjust = 1)
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.title = element_blank()) +
        ylab(paste(analy_full, unit, sep = ""))
      if (faceted) plot <- plot + facet_wrap(~partition)
      
      print(plot)
      dev.off()
    }
    
    # with error bars, with comb_df and working dodging
    pd <- position_dodge(width = 0.2)
    
    if (dataset == "clinical") {
      # setting up background rects
      data_min <- min(comb_df$mean_parts - comb_df$std_error) # comb with error makes sure rects cover error bars too
      data_max <- max(comb_df$mean_parts + comb_df$std_error)
      # data_min <- data_mins[[analy]]
      # data_max <- data_maxs[[analy]]
      if (data_min >= clin_maxs[analy] || data_max <= clin_mins[analy]) { # case of no normal_range
        ranges <- data.frame(ystart = data_min,
                             yend = data_max,
                             col = "outside_normal")
      } else if (data_min >= clin_mins[analy] && data_max <= clin_maxs[analy]) { # case of normal_range only
        ranges <- data.frame(ystart = data_min,
                             yend = data_max,
                             col = "normal_range")
      } else if (data_min >= clin_mins[analy]) { # case of normal_range + outside_normal upper
        ranges <- data.frame(ystart = c(data_min, clin_maxs[analy]),
                             yend = c(clin_maxs[analy], data_max),
                             col = c("normal_range", "outside_normal"))
      } else if (data_max <= clin_maxs[analy]) { # case of normal_range + outside_normal lower
        ranges <- data.frame(ystart = c(data_min, clin_mins[analy]),
                             yend = c(clin_mins[analy], data_max),
                             col = c("outside_normal", "normal_range"))
      } else { # case normal_range + outside_normal on both sides
        ranges <- data.frame(ystart = c(data_min, clin_mins[analy], clin_maxs[analy]),
                             yend = c(clin_mins[analy], clin_maxs[analy], data_max),
                             col = c("outside_normal", "normal_range", "outside_normal"))
      }
    }
    
    pdf(file.path(graph_dir, paste(analy, "_", dataset, "_Comb_Avg.pdf", sep = "")), width = 9, height = 6)
    
    # normalizing with each fiber baseline
    # assumed sorted, so baseline is first timepoint of each analyte-fiber set, and fibers adjacent
    
    plot <- comb_df %>% ggplot()
    if (dataset == "clinical") plot <- plot +
      geom_rect(data = ranges, aes(ymin = ystart, ymax = yend, xmin = -Inf, xmax = Inf, fill = col), alpha = 0.4) +
      scale_fill_manual(values = c(outside_normal = "#ffcccb", normal_range = "#c6ff95"))
    plot <- plot +  
      geom_line(position = pd, aes(x = point, y = mean_parts, group = fiber, color = fiber)) +
      geom_errorbar(aes(x = point, ymin = mean_parts - std_error, ymax = mean_parts + std_error,
                        group = fiber, color = fiber),
                    width = .5, position = pd) +
      theme(
        axis.text.x = element_text(size=10, angle=45, hjust = 1)) +
        ylab(paste(analy_full, unit, sep = ""))
    
    print(plot)
    dev.off()
  }
}