# investigates ONE analyte in an omic dataset and makes a fiber-comb graph and individual graphs
#   per fiber of the participants' responses
# adds overlay of fiber average on partitioned data
# includes normal ranges, units, full names for clinicals
# includes full names for proteomics
# "cleanup" changes for faceted:
#   removed x-axis tick labels, x-axis title ("point"), y-axis "baselines normalized to mean," legend, different id colors 

# @dataset      e.g. genef, pcl, metaphlan, cytokine, clinical, lipids, proteomics, metabolomics
# @selected     analyte(s), default = all analytes
# @norm         whether normalized with all baselines @ 0 (precomputed)

# @comb_only              only print combined-fiber graphs
# @faceted                faceted instead of overlayed ids
# @filled                 "fill in" missing baselines with avg of present baselines (precomputed)
# @partition              how to divide the data, default = ids separate
# @responder_label        artificially create binary division, which partition label to be treated as responder
# @partition_method       partitioning methods, current options: transparency, aggregate
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

analyte_investigate_overlay <- function(dataset, selected = all_analytes, norm = T, 
                                        comb_only = F, faceted = F, filled = F, partition = NA, responder_label = "responder",
                                        partition_method = "transparency",
                                        only = ids, without = c(), omit_x_axis = F, omit_units = F,
                                        fibers = all_fibers, file_prefix = "Tidy_Full", graph_dir = NA,
                                        responder_partition = T, label_responders = T, desc = "") {
  require(tidyverse)
  require(plotrix)
  
  if (filled) file_prefix <- paste(file_prefix, "Filled", sep = "_")

  all_fibers = scan(file.path("Metadata", "Fibers.tsv"), character(), quote = '', sep = "\t", quiet = T)
  ids = scan(file.path("Metadata", "Ids.tsv"), character(), quote = '', sep = "\t", quiet = T)
  load(file.path("Data", file_prefix, paste(file_prefix, fibers[1], dataset, "df.RData", sep = "_")))
  all_analytes <- tidy_df$analyte %>% unique()
  
  if (is.na(graph_dir)) {
    dir_source <- "Graphs"
    dir_name <- "Analyte"
    if (norm) dir_name <- paste(dir_name, "Norm", sep = "_")
    if (comb_only) dir_name <- paste(dir_name, "Comb_Only", sep = "_")
    if (faceted) dir_name <- paste(dir_name, "Faceted", sep = "_")
    if (filled) dir_name <- paste(dir_name, "Filled", sep = "_")
    if (!is.na(partition)) dir_name <- paste(dir_name, "Partition_Overlay", sep = "_")
    if (length(selected) == 1) dir_name <- paste(dir_name, selected, sep = "_")
    if (length(fibers) == 1) dir_name <- paste(dir_name, fibers, sep = "_")
    dir_name <- paste(dir_name, dataset, sep = "_")
    if (desc != "") dir_name <- paste(dir_name, desc, sep = "_")
    graph_dir <- file.path(dir_source, dir_name)
  }
  if (!dir.exists(graph_dir))
    dir.create(graph_dir)
  
  if (dataset == "clinical") {
    # add units
    clin_names <- scan(file.path("Metadata", "Clin_Names_Safe.tsv"), character(), quote = '', sep = "\t", quiet = T)
    clin_names_full <- scan(file.path("Metadata", "Clin_Names.tsv"), character(), quote = '', sep = "\t", quiet = T)
    clin_units <- scan(file.path("Metadata", "Clin_Units.tsv"), character(), quote = '', sep = "\t", quiet = T)
    names(clin_units) <- clin_names
    names(clin_names_full) <- clin_names
    
    # also reading for normal_range ranges
    clin_mins <- scan(file.path("Metadata", "Clin_Healthy_Ranges_min.tsv"), double(), quote = '', sep = "\t", quiet = T)
    clin_maxs <- scan(file.path("Metadata", "Clin_Healthy_Ranges_max.tsv"), double(), quote = '', sep = "\t", quiet = T)
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
    names(partition_vec) <- partition_df$id
  }
  
  # if no partition, responder partitioning is impossible
  if (is.na(partition)) responder_partition = F
  # transparency only if chosen and responder_partition available and not faceted
  if (responder_partition && !faceted && partition_method == "transparency") {
    transparency = T
  } else transparency = F
  # aggregate only if chosen and responder_partition available and not faceted
  if (responder_partition && !faceted && partition_method == "aggregate") {
    aggregate = T
  } else aggregate = F
  # if no responder partitioning or using aggregate, do not label responders
  if (!responder_partition || aggregate) label_responders = F
  
  # dodge for aggreg or combined fiber graphs
  pd <- position_dodge(width = 0.2)
  
  for (analy in selected) {
    
    comb_df <- data.frame()
  
    # setting units
    if (omit_units) {
      unit <- ""
    } else if (dataset == "clinical") {
      if (is.na(clin_units[analy])) {
        unit <- ""
      } else {
        unit <- paste(" (", clin_units[analy], ")", sep = "")
      }
      # if (norm) {
      #   unit <- paste(unit, ",\nbaselines normalized to mean", sep = "")
      # }
    } else {
      if (norm) {
        unit <- "\n(log2 normalized expression)"
      } else {
        unit <- "\n(log2 expression)"
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
      load(file.path("Data", file_prefix, paste(file_prefix, fiber, dataset, "df.RData", sep = "_")))
      tidy_df <- tidy_df %>% filter(analyte == analy) %>% filter(id %in% only) %>% filter(!(id %in% without))
      
      # choosing norm or not
      if (norm) {
        tidy_df <- na.omit(tidy_df) # for missing baselines creating NAs in normalization
        tidy_df <- tidy_df %>% mutate(value = renorm_val)
      } else {
        tidy_df <- tidy_df %>% mutate(value = val)
      }
      
      # graphing ids together, averaged
      tmp_df <- tidy_df %>%
        group_by(point) %>%
        summarize(mean_parts = mean(value), std_error = std.error(value)) %>%
        ungroup()

      # saving dfs
      fiber_av_df <- tmp_df
      tmp_df <- cbind(tmp_df, "fiber" = rep(fiber, nrow(tmp_df)))
      comb_df <- rbind(comb_df, tmp_df)
      
      # check if we want individual (non-comb) graphs
      if (comb_only) next
      
      if (dataset == "clinical") {
        # setting up background rects
        data_min <- min(tidy_df$value)
        data_max <- max(tidy_df$value)
        
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
        tidy_df$partition = tidy_df$id
      } else {
        tidy_df$partition = partition_vec[tidy_df$id]
      }
      
      # responder partitioning: find responders
      if (responder_partition) {
        responders = tidy_df[which(partition_vec[tidy_df$id] == responder_label), "id"] %>% unique()
      }
      
      # transparency
      if (transparency) {
        tidy_df$alph <- partition_vec[tidy_df$id] == responder_label
        # generate random colors for responders
        # colors_per_id = rainbow(length(ids))[sample(1:length(ids), length(ids))]
        colors_per_id = seq(15, 375, length = length(ids))
        colors_per_id = hcl(h = colors_per_id, l = 65, c = 100)[1:length(ids)]
        names(colors_per_id) <- c(ids)
        colors_per_id[partition_vec[ids] != responder_label] = "#C2C5CC"
        
        responders <- ids[partition_vec[ids] == responder_label]
      } else if (aggregate) {
        tidy_df$partition <- partition_vec[tidy_df$id]
        tidy_df <- tidy_df %>%
          group_by(point, partition) %>%
          summarize(mean_partition = mean(value), std_error = std.error(value)) %>% # std error cannot be based on norm
          ungroup()
      }
      
      # for grouping
      fiber_av_df$id <- "average"
      
      # graphing ids separately
      pdf(file.path(graph_dir, paste(analy, dataset, fiber, "Ids.pdf", sep = "_")), width = 6, height = 4)
      
      plot <- tidy_df %>% ggplot()
      if (dataset == "clinical") plot <- plot +
        geom_rect(data = ranges, aes(ymin = ystart, ymax = yend, xmin = -Inf, xmax = Inf, fill = col), alpha = 0.4) +
        scale_fill_manual(name = "range", values = c(outside_normal = "#ffcccb", normal_range = "#c6ff95"))
      # manage colors
      if (faceted && is.na(partition)) plot <- plot +  # faceted and no partition: no colors
        geom_line(aes(x = point, y = value, group = id))
      else { # otherwise, colors
        if (transparency) plot <- plot +
          geom_line(aes(x = point, y = value, group = id, alpha = alph, color = id))
        else if (aggregate) plot <- plot +
          geom_line(aes(x = point, y = mean_partition, group = partition, color = partition), position = pd) +
          geom_errorbar(aes(x = point, ymin = mean_partition - std_error, ymax = mean_partition + std_error,
                            group = partition, color = partition),
                        width = .5, position = pd)
        else plot <- plot +
          geom_line(aes(x = point, y = value, group = id, color = id))
      }
      # theme management
      if (omit_x_axis) plot <- plot +
        theme(
          # axis.text.x = element_text(size=10, angle=45, hjust = 1)
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none") # remove both id legend and potential alpha legend
      else plot <- plot +
        theme(
          axis.text.x = element_text(size=10, angle=45, hjust = 1)) +
        guides(alpha = F) # remove only potential alpha legend
      plot <- plot +
        ylab(paste(analy_full, unit, sep = ""))
      # optionals
      if (transparency) {
        plot <- plot + scale_alpha_manual(values = c(.3, 1), labels = NULL) +
          scale_color_manual(values = c(colors_per_id, "average" = "black"), labels = c(responders, "fiber_average"), breaks = c(responders, "average"), name = responder_label)
      }
      # if (label_responders) plot <- plot + scale_color_discrete(breaks = responders, name = responder_label)
      if (faceted) plot <- plot + facet_wrap(~partition)
      
      # add overlay
      # print(c(colors_per_id, average = "black"))
      # print(names(c(colors_per_id, average = "black")))
      # print(c(responders, "average"))
      plot <- plot +
        geom_line(data = fiber_av_df, position = pd, aes(x = point, y = mean_parts, group = 1, color = "average")) +
        geom_errorbar(data = fiber_av_df, aes(x = point, ymin = mean_parts - std_error, ymax = mean_parts + std_error, 
                                              group = 1, color = "average"), width = .5, position = pd)
      
      print(plot)
      dev.off()
    }
    
    if (length(fibers) == 1) return("Only 1 fiber, so no comb graph")
    
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
    
    pdf(file.path(graph_dir, paste(analy, "_", dataset, "_Comb_Avg.pdf", sep = "")), width = 6, height = 4)
    
    # normalizing with each fiber baseline
    # assumed sorted, so baseline is first timepoint of each analyte-fiber set, and fibers adjacent
    
    plot <- comb_df %>% ggplot()
    if (dataset == "clinical") plot <- plot +
      geom_rect(data = ranges, aes(ymin = ystart, ymax = yend, xmin = -Inf, xmax = Inf, fill = col), alpha = 0.4) +
      scale_fill_manual(name = "range", values = c(outside_normal = "#ffcccb", normal_range = "#c6ff95"))
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
