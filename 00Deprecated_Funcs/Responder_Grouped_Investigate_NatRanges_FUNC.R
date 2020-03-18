# investigate ONE analy as 3 groups on 1 graph, with top, middle, bottom 6 responders
# responder val calculated with mean of 10, 20, 30
# @precondition there are 18 ids

# @dataset    e.g. genef, pcl, metaphlan, cytokine, clinical, lipids, proteomics, metabolomics
# @analy      desired analyte in the dataset
# @norm       whether normalized with all baselines @ 0

responder_grouped_investigate <- function(dataset, analy, norm) {

  require(tidyverse)
  require(plotrix)
  
  fibers = scan(file.path("Metadata", "Fibers.tsv"), character(), quote = '', sep = "\t", quiet = T)
  file_prefix <- "Tidy_Full"
  
  # dirs
  if (norm) {
    dir_name <- "Responder_Grouped_Graphs_Norm"
  } else {
    dir_name <- "Responder_Grouped_Graphs"
  }
  graph_dir <- paste(dir_name, analy, dataset, sep = "_")
  if (!dir.exists(graph_dir)) {
    dir.create(graph_dir)
  }
  
  # unit setup
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
  
  # setting units
  # setting units
  if (dataset == "clinical") {
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
    # getting data for grouping, updated: uses natural norm (Tidy_Full)
    # THUS ONLY FOR CLINICALS
    # load(paste(file_prefix, fiber, dataset, "median_df.RData", sep = "_"))
    load(paste("Tidy_Full", fiber, dataset, "df.RData", sep = "_"))
    tidy_cyto2 <- tidy_df
    
    # chol_df <- tidy_cyto2 %>% filter(point == "30" | point == "Baseline") %>%
    #   select(id, point, analy) %>%
    #   spread(point, analy, sep = "_") %>%
    #   mutate(diff = point_30 - point_Baseline) %>%
    #   select(id, diff) %>%
    #   arrange(diff)
    # below ver is for median df, ungathered
    # resp_df <- tidy_cyto2 %>% filter(point == "30" | point == "20" | point == "10") %>%
    #   select(id, point, analy) %>%
    #   spread(point, analy, sep = "_") %>%
    #   mutate(avg = rowMeans(dplyr::select(., point_30, point_20, point_10), na.rm = T)) %>%
    #   select(id, avg) %>%
    #   arrange(avg)
    if (norm) {
      resp_df <- tidy_cyto2 %>% filter(point == "30" | point == "20" | point == "10") %>%
        filter(analyte == analy) %>%
        select(-c(val, analyte)) %>% # necessary for clean/proper spread
        spread(point, norm_val, sep = "_") %>%
        mutate(avg = rowMeans(dplyr::select(., point_30, point_20, point_10), na.rm = T)) %>%
        select(id, avg) %>%
        arrange(avg)
    } else {
      resp_df <- tidy_cyto2 %>% filter(point == "30" | point == "20" | point == "10") %>%
        filter(analyte == analy) %>%
        select(-c(norm_val, analyte)) %>%
        spread(point, val, sep = "_") %>%
        mutate(avg = rowMeans(dplyr::select(., point_30, point_20, point_10), na.rm = T)) %>%
        select(id, avg) %>%
        arrange(avg)
    }
    
    # making groups
    response <- factor(levels = c("top", "middle", "bottom"))
    resp_df$resp <- levels(response)[(0:(nrow(resp_df) - 1) %/% 6) + 1] # splitting into 3 groups
    resp_vec <- resp_df$resp
    names(resp_vec) <- resp_df$id # vector more convenient for checking
    
    # getting data for plot
    # load(paste(file_prefix, "Log", fiber, dataset, "df2.RData", sep = "_")) # flaw: currently the gathered (df2)
    #                                                                  # datasets are all named log, even clin
    tidy_df <- tidy_cyto2 %>% mutate(resp = resp_vec[id])
    if (norm) {
      tidy_df <- na.omit(tidy_df) # for missing baselines creating NAs in normalization
    }
    tmp_df <- tidy_df %>% filter(analyte == analy)
    
    # averaging across resp groups
    if (norm) {
      tmp_df <- tmp_df %>%
        group_by(point, resp) %>%
        dplyr::summarise(mean_parts = mean(norm_val), std_error = std.error(norm_val)) %>% # std error cannot be based on norm
        ungroup()
    } else {
      tmp_df <- tmp_df %>%
        group_by(point, resp) %>%
        dplyr::summarise(mean_parts = mean(val), std_error = std.error(val)) %>% # std error cannot be based on norm
        ungroup()
    }
    
    # settingst up ranges
    if (dataset == "clinical") {
      # setting up background rects
      data_min <- min(tmp_df$mean_parts - tmp_df$std_error) # comb with error makes sure rects cover error bars too
      data_max <- max(tmp_df$mean_parts + tmp_df$std_error)
      
      # data_min <- data_mins[[analy]]
      # data_max <- data_maxs[[analy]]
      if (data_min >= clin_maxs[analy] || data_max <= clin_mins[analy]) { # case of no healthy
        ranges <- data.frame(ystart = data_min,
                             yend = data_max,
                             col = "outside_normal")
      } else if (data_min >= clin_mins[analy] && data_max <= clin_maxs[analy]) { # case of healthy only
        ranges <- data.frame(ystart = data_min,
                             yend = data_max,
                             col = "normal_range")
      } else if (data_min >= clin_mins[analy]) { # case of healthy + unhealthy upper
        ranges <- data.frame(ystart = c(data_min, clin_maxs[analy]),
                             yend = c(clin_maxs[analy], data_max),
                             col = c("normal_range", "outside_normal"))
      } else if (data_max <= clin_maxs[analy]) { # case of healthy + unhealthy lower
        ranges <- data.frame(ystart = c(data_min, clin_mins[analy]),
                             yend = c(clin_mins[analy], data_max),
                             col = c("outside_normal", "normal_range"))
      } else { # case healthy + unhealthy on both sides
        ranges <- data.frame(ystart = c(data_min, clin_mins[analy], clin_maxs[analy]),
                             yend = c(clin_mins[analy], clin_maxs[analy], data_max),
                             col = c("outside_normal", "normal_range", "outside_normal"))
      }
    }
    
    # with error bars, use dodging
    pd <- position_dodge(width = 0.2)
    
    # plotting
    pdf(file.path(graph_dir, paste(analy, dataset, fiber, "Response_Grouped_Graph.pdf", sep = "_")), width = 6, height = 4)
    
    if (dataset == "clinical") {
      print(tmp_df %>%
              ggplot(position = pd) +
              geom_rect(data = ranges, aes(ymin = ystart, ymax = yend, xmin = -Inf, xmax = Inf, fill = col), alpha = 0.4) +
              scale_fill_manual(values = c(outside_normal = "#ffcccb", normal_range = "#c6ff95")) +
              geom_line(position = pd, aes(x = point, y = mean_parts, group = resp, color = resp)) +
              geom_errorbar(aes(x = point, ymin = mean_parts - std_error, ymax = mean_parts + std_error,
                                group = resp, color = resp),
                            width = .5, position = pd) +
              theme(
                axis.text.x = element_text(size=10, angle=45, hjust = 1)) +
              ylab(paste(analy_full, unit, sep = "")))
    } else {
      print(tmp_df %>%
              ggplot(position = pd, aes(x = point, y = mean_parts, group = resp, color = resp)) +
              geom_line(position = pd) +
              geom_errorbar(aes(x = point, ymin = mean_parts - std_error, ymax = mean_parts + std_error),
                            width = .5, position = pd) +
              theme(
                axis.text.x = element_text(size=10, angle=45, hjust = 1)) +
              ylab(paste(analy_full, unit, sep = "")))
    }
    dev.off()
  }
}