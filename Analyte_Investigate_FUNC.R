# investigates ONE analyte in an omic/dataset and makes a fiber-comb graph and individual graphs
#   per fiber showing the participants, so 4 graphs total
# renorms centered at part means
# adds normal_range ranges for clinicals
# updated with units, takes care of protein full names and clinical units in particular

# @dataset    choices here are genef, pcl, metaphlan, cytokine, clinical, lipids, proteomics, 
#             metabolomics
# @selected   vector of analyte string(s), default = all analytes
# @renorm     whether renormalized with all baselines @ 0
# @list_name  optional, for documentation, needs to describe clearly

# "quick and dirty" way to get all analytes, assume arabinoxylan file
load(paste("Tidy_Full", "Arabinoxylan", dataset, "df.RData", sep = "_"))
all_analytes <- tidy_df$analyte %>% unique()

analyte_investigate_list_natural_ranges <- function(dataset, selected = all_analytes, renorm, list_name = "") {

  require(tidyverse)
  require(plotrix)
  
  fibers = c("Arabinoxylan","LCInulin","Mix")
  
  # order matters as this is used for factor levels
  timepoints <- c("Baseline", "10", "20", "30",
                  "WashoutD3", "WashoutD10", "WashoutFinal")
  
  if (renorm) {
    dir_name <- "Analyte_Analyses_Natural_Ranges_Renorm"
  } else {
    dir_name <- "Analyte_Analyses_Natural_Ranges"
  }
  
  graph_dir <- paste(dir_name, analy, dataset, list_name, sep = "_")
  if (!dir.exists(graph_dir))
    dir.create(graph_dir)
  
  if (dataset == "clinical") {
    # add units
    clin_names <- scan("Clin_Names_Safe.tsv", character(), quote = '', sep = "\t")
    clin_names_full <- scan("Clin_Names.tsv", character(), quote = '', sep = "\t")
    clin_units <- scan("Clin_Units.tsv", character(), quote = '', sep = "\t")
    names(clin_units) <- clin_names
    names(clin_names_full) <- clin_names
    
    # also reading for normal_range ranges
    clin_mins <- scan("Clin_Healthy_Ranges_min.tsv", double(), quote = '', sep = "\t")
    clin_maxs <- scan("Clin_Healthy_Ranges_max.tsv", double(), quote = '', sep = "\t")
    names(clin_mins) <- clin_names
    names(clin_maxs) <- clin_names
  }
  
  if (dataset == "proteomics") {
    # loading data
    load("Proteomics.RData")
    proteo <- Proteomics_data_1pept_Log2_filt80perc_ImptND_Combat
    # class(proteo$Protein.name) <- "character"
    proteo$Protein.name <- as.character(proteo$Protein.name)
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
      if (renorm) {
        unit <- paste(unit, ",\nbaselines normalized to mean", sep = "")
      }
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
    anl <- analy
    
    for (fiber in fibers) {
      print(fiber)
      
      # loading data
      # load(paste("Tidy_Normalized_Log", fiber, dataset, "df2.RData", sep = "_")
      load(paste("Tidy_Full", fiber, dataset, "df.RData", sep = "_"))
      tidy_cyto <- tidy_df %>% filter(analyte == analy)
      
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
        # data_min <- data_mins[[anl]]
        # data_max <- data_maxs[[anl]]
        if (data_min >= clin_maxs[anl] || data_max <= clin_mins[anl]) { # case of no normal_range
          ranges <- data.frame(ystart = data_min,
                               yend = data_max,
                               col = "outside_normal")
        } else if (data_min >= clin_mins[anl] && data_max <= clin_maxs[anl]) { # case of normal_range only
          ranges <- data.frame(ystart = data_min,
                               yend = data_max,
                               col = "normal_range")
        } else if (data_min >= clin_mins[anl]) { # case of normal_range + outside_normal upper
          ranges <- data.frame(ystart = c(data_min, clin_maxs[anl]),
                               yend = c(clin_maxs[anl], data_max),
                               col = c("normal_range", "outside_normal"))
        } else if (data_max <= clin_maxs[anl]) { # case of normal_range + outside_normal lower
          ranges <- data.frame(ystart = c(data_min, clin_mins[anl]),
                               yend = c(clin_mins[anl], data_max),
                               col = c("outside_normal", "normal_range"))
        } else { # case normal_range + outside_normal on both sides
          ranges <- data.frame(ystart = c(data_min, clin_mins[anl], clin_maxs[anl]),
                               yend = c(clin_mins[anl], clin_maxs[anl], data_max),
                               col = c("outside_normal", "normal_range", "outside_normal"))
        }
      }
      
      # graphing ids separately
      if (dataset == "clinical") {
        if (renorm) {
          pdf(file.path(graph_dir, paste(analy, dataset, fiber, "Ids.pdf", sep = "_")), width = 9, height = 6)
          print(tidy_cyto %>%
                  ggplot() +
                  geom_rect(data = ranges, aes(ymin = ystart, ymax = yend, xmin = -Inf, xmax = Inf, fill = col), alpha = 0.4) +
                  scale_fill_manual(values = c(outside_normal = "#ffcccb", normal_range = "#c6ff95")) +
                  # geom_point() +
                  geom_line(aes(x = point, y = renorm_val, group = id, color = id)) +
                  theme(
                    axis.text.x = element_text(size=10, angle=45, hjust = 1)) +
                  ylab(paste(analy_full, unit, sep = "")))
          dev.off()
        } else {
          pdf(file.path(graph_dir, paste(analy, dataset, fiber, "Ids.pdf", sep = "_")), width = 9, height = 6)
          print(tidy_cyto %>%
                  ggplot() +
                  geom_rect(data = ranges, aes(ymin = ystart, ymax = yend, xmin = -Inf, xmax = Inf, fill = col), alpha = 0.4) +
                  scale_fill_manual(values = c(outside_normal = "#ffcccb", normal_range = "#c6ff95")) +
                  # geom_point() +
                  geom_line(aes(x = point, y = val, group = id, color = id)) +
                  theme(
                    axis.text.x = element_text(size=10, angle=45, hjust = 1)) +
                  ylab(paste(analy_full, unit, sep = "")))
          dev.off()
        }
      } else {
        if (renorm) {
          pdf(file.path(graph_dir, paste(analy, dataset, fiber, "Ids.pdf", sep = "_")), width = 9, height = 6)
          print(tidy_cyto %>%
                  ggplot() +
                  # geom_rect(data = ranges, aes(ymin = ystart, ymax = yend, xmin = -Inf, xmax = Inf, fill = col), alpha = 0.4) +
                  # scale_fill_manual(values = c(outside_normal = "#ffcccb", normal_range = "#c6ff95")) +
                  # geom_point() +
                  geom_line(aes(x = point, y = renorm_val, group = id, color = id)) +
                  theme(
                    axis.text.x = element_text(size=10, angle=45, hjust = 1)) +
                  ylab(paste(analy_full, unit, sep = "")))
          dev.off()
        } else {
          pdf(file.path(graph_dir, paste(analy, dataset, fiber, "Ids.pdf", sep = "_")), width = 9, height = 6)
          print(tidy_cyto %>%
                  ggplot() +
                  # geom_rect(data = ranges, aes(ymin = ystart, ymax = yend, xmin = -Inf, xmax = Inf, fill = col), alpha = 0.4) +
                  # scale_fill_manual(values = c(outside_normal = "#ffcccb", normal_range = "#c6ff95")) +
                  # geom_point() +
                  geom_line(aes(x = point, y = val, group = id, color = id)) +
                  theme(
                    axis.text.x = element_text(size=10, angle=45, hjust = 1)) +
                  ylab(paste(analy_full, unit, sep = "")))
          dev.off()
        }
      }
    }
    
    # with error bars, with comb_df and working dodging
    pd <- position_dodge(width = 0.2)
    
    
    if (dataset == "clinical") {
      # setting up background rects
      data_min <- min(comb_df$mean_parts - comb_df$std_error) # comb with error makes sure rects cover error bars too
      data_max <- max(comb_df$mean_parts + comb_df$std_error)
      # data_min <- data_mins[[anl]]
      # data_max <- data_maxs[[anl]]
      if (data_min >= clin_maxs[anl] || data_max <= clin_mins[anl]) { # case of no normal_range
        ranges <- data.frame(ystart = data_min,
                             yend = data_max,
                             col = "outside_normal")
      } else if (data_min >= clin_mins[anl] && data_max <= clin_maxs[anl]) { # case of normal_range only
        ranges <- data.frame(ystart = data_min,
                             yend = data_max,
                             col = "normal_range")
      } else if (data_min >= clin_mins[anl]) { # case of normal_range + outside_normal upper
        ranges <- data.frame(ystart = c(data_min, clin_maxs[anl]),
                             yend = c(clin_maxs[anl], data_max),
                             col = c("normal_range", "outside_normal"))
      } else if (data_max <= clin_maxs[anl]) { # case of normal_range + outside_normal lower
        ranges <- data.frame(ystart = c(data_min, clin_mins[anl]),
                             yend = c(clin_mins[anl], data_max),
                             col = c("outside_normal", "normal_range"))
      } else { # case normal_range + outside_normal on both sides
        ranges <- data.frame(ystart = c(data_min, clin_mins[anl], clin_maxs[anl]),
                             yend = c(clin_mins[anl], clin_maxs[anl], data_max),
                             col = c("outside_normal", "normal_range", "outside_normal"))
      }
    }
    
    pdf(file.path(graph_dir, paste(analy, "_", dataset, "Comb_Avg.pdf", sep = "")), width = 9, height = 6)
    
    # normalizing with each fiber baseline
    # assumed sorted, so baseline is first timepoint of each analyte-fiber set, and fibers adjacent
    
    if (dataset == "clinical") {
      print(comb_df %>%
              ggplot() + #ggplot(position = pd) +
              geom_rect(data = ranges, aes(ymin = ystart, ymax = yend, xmin = -Inf, xmax = Inf, fill = col), alpha = 0.4) +
              scale_fill_manual(values = c(outside_normal = "#ffcccb", normal_range = "#c6ff95")) +
              geom_line(position = pd, aes(x = point, y = mean_parts, group = fiber, color = fiber)) +
              geom_errorbar(aes(x = point, ymin = mean_parts - std_error, ymax = mean_parts + std_error,
                                group = fiber, color = fiber),
                            width = .5, position = pd) +
              theme(
                axis.text.x = element_text(size=10, angle=45, hjust = 1)) +
              ylab(paste(analy_full, unit, sep = "")))
      dev.off()
    } else {
      print(comb_df %>%
              ggplot() + #ggplot(position = pd) +
              # geom_rect(data = ranges, aes(ymin = ystart, ymax = yend, xmin = -Inf, xmax = Inf, fill = col), alpha = 0.4) +
              # scale_fill_manual(values = c(outside_normal = "#ffcccb", normal_range = "#c6ff95")) +
              geom_line(position = pd, aes(x = point, y = mean_parts, group = fiber, color = fiber)) +
              geom_errorbar(aes(x = point, ymin = mean_parts - std_error, ymax = mean_parts + std_error,
                                group = fiber, color = fiber),
                            width = .5, position = pd) +
              theme(
                axis.text.x = element_text(size=10, angle=45, hjust = 1)) +
              ylab(paste(analy_full, unit, sep = "")))
      dev.off()
    }
  }
}