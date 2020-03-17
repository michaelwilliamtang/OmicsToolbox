# investigates ONE analyte in an omic/dataset and makes a fiber-comb graph and individual graphs
#   per fiber showing the participants, so 4 graphs total
# renorms centered at part means
# adds healthy ranges for clinicals
# updated with units, takes care of protein full names and clinical units in particular

# @dataset    choices here are genef, pcl, metaphlan, cytokine, clinical, lipids, proteomics, 
#             metabolomics
# @analy      string of which analyte in the dataset you want
# @renorm     whether renormalized with all baselines @ 0

# sex split: (g/dL)
# m: 13.5-17.5
# f: 12-15.5
# written as function, but ranges provided are specifically for hemoglobin
# if generalizable need arises, will compile data file of similar sex ranges, for now hard-coded to hemoglobin
# THIS VER: updated with sex-combined Avg graph within fibers (compare to old fiber-combined Avg graph)

analyte_investigate_natural_ranges_sex <- function(dataset, analy, renorm) {

  require(tidyverse)
  require(plotrix)
  
  # sex data on participants
  part_sex <- read.csv("part_sex.csv", stringsAsFactors = F)
  males <- part_sex$participant[which(part_sex$Sex == 'M')]
  females <- part_sex$participant[which(part_sex$Sex == 'F')]
  
  # specific data on hemo and other params
  # dataset <- "clinical"
  # analy <- "Hemoglobin"
  # renorm <- T
  sex_ranges <- list()
  sex_ranges[["male"]] <- c(13.5, 17.5)
  sex_ranges[["female"]] <- c(12, 15.5)
  sexes <- c("male", "female")
  
  fibers = c("Arabinoxylan","LCInulin","Mix")
  
  # order matters as this is used for factor levels
  timepoints <- c("Baseline", "10", "20", "30",
                  "WashoutD3", "WashoutD10", "WashoutFinal")
  
  if (renorm) {
    dir_name <- "Analyte_Analyses_Natural_Ranges_Sex_Renorm"
  } else {
    dir_name <- "Analyte_Analyses_Natural_Ranges_Sex"
  }
  
  graph_dir <- paste(dir_name, analy, dataset, sep = "_")
  if (!dir.exists(graph_dir))
    dir.create(graph_dir)
  
  # comb_df <- data.frame()
  # comb_sexes <- list()
  # comb_sexes[["male"]] <- data.frame()
  # comb_sexes[["female"]] <- data.frame()
  
  if (dataset == "clinical") {
    # add units
    clin_names <- scan("Clin_Names_Safe.tsv", character(), quote = '', sep = "\t")
    clin_names_full <- scan("Clin_Names.tsv", character(), quote = '', sep = "\t")
    clin_units <- scan("Clin_Units.tsv", character(), quote = '', sep = "\t")
    names(clin_units) <- clin_names
    names(clin_names_full) <- clin_names
    
    # also reading for healthy ranges
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
    unit <- " (log2 normalized expression)"
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

    comb_sex <- data.frame()
    
    # loading data
    # load(paste("Tidy_Normalized_Log", fiber, dataset, "df2.RData", sep = "_")
    load(paste("Tidy_Full", fiber, dataset, "df.RData", sep = "_"))
    tidy_cyto <- tidy_df %>% filter(analyte == analy)
    
    # split by sex
    tidy_sex <- list()
    tidy_sex[["male"]] <- tidy_cyto %>% filter(id %in% males)
    tidy_sex[["female"]] <- tidy_cyto %>% filter(id %in% females)
    
    ##### get participant means by sex and add to comb dfs
    for (sx in sexes) {
      tidy_cyto <- tidy_sex[[sx]]
  
      # graphing ids together, averaged
      if (renorm) {
        # SPECIAL: create new custom renorms within each sex
        source("Renormalize_Natural2.R")
        tidy_cyto <- renormalize_natural(tidy_cyto %>% select(-renorm_val)) # select anti is optional, included as anti-bug guarantee for renorm_val overwrite
  
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
    
      # saving dfs for avg graph
      tmp_cyto <- cbind(tmp_cyto, "sex" = rep(sx, nrow(tmp_cyto)))
      comb_sex <- rbind(comb_sex, tmp_cyto)
  
      ##### setting up ranges based on sex
      # tidy_cyto <- tidy_sex[[sx]]
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
        if (data_min >= sex_ranges[[sx]][2] || data_max <= sex_ranges[[sx]][1]) { # case of no normal_range
          ranges <- data.frame(ystart = data_min,
                               yend = data_max,
                               col = "outside_normal")
        } else if (data_min >= sex_ranges[[sx]][1] && data_max <= sex_ranges[[sx]][2]) { # case of normal_range only
          ranges <- data.frame(ystart = data_min,
                               yend = data_max,
                               col = "normal_range")
        } else if (data_min >= sex_ranges[[sx]][1]) { # case of normal_range + outside_normal upper
          ranges <- data.frame(ystart = c(data_min, sex_ranges[[sx]][2]),
                               yend = c(sex_ranges[[sx]][2], data_max),
                               col = c("normal_range", "outside_normal"))
        } else if (data_max <= sex_ranges[[sx]][2]) { # case of normal_range + outside_normal lower
          ranges <- data.frame(ystart = c(data_min, sex_ranges[[sx]][1]),
                               yend = c(sex_ranges[[sx]][1], data_max),
                               col = c("outside_normal", "normal_range"))
        } else { # case normal_range + outside_normal on both sides
          ranges <- data.frame(ystart = c(data_min, sex_ranges[[sx]][1], sex_ranges[[sx]][2]),
                               yend = c(sex_ranges[[sx]][1], sex_ranges[[sx]][2], data_max),
                               col = c("outside_normal", "normal_range", "outside_normal"))
        }
      }
      # graphing ids separately
      if (renorm) {
        pdf(file.path(graph_dir, paste(analy, dataset, fiber, sx, "Ids.pdf", sep = "_")), width = 9, height = 6)
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
        pdf(file.path(graph_dir, paste(analy, dataset, fiber, sx, "Ids.pdf", sep = "_")), width = 9, height = 6)
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
    }
    
    # avg graph between sexes, omit ranges since they are different between sexes
    pdf(file.path(graph_dir, paste(analy, dataset, fiber, "Avg.pdf", sep = "_")), width = 9, height = 6)
    pd <- position_dodge(width = 0.2)

    print(comb_sex %>%
        ggplot() + #ggplot(position = pd) +
        # geom_rect(data = ranges, aes(ymin = ystart, ymax = yend, xmin = -Inf, xmax = Inf, fill = col), alpha = 0.4) +
        # scale_fill_manual(values = c(outside_normal = "#ffcccb", normal_range = "#c6ff95")) +
        geom_line(position = pd, aes(x = point, y = mean_parts, group = sex, color = sex)) +
        geom_errorbar(aes(x = point, ymin = mean_parts - std_error, ymax = mean_parts + std_error,
                          group = sex, color = sex),
                      width = .5, position = pd) +
        theme(
          axis.text.x = element_text(size=10, angle=45, hjust = 1)) +
        ylab(paste(analy_full, unit, sep = "")))
    dev.off()
  }
  ###############################
  # NOTE: omit ranges since they are different between sexes
  
  # with error bars, with comb_df and working dodging
  # pd <- position_dodge(width = 0.2)
  
  # for (sx in sexes) {
  #   comb_df <- comb_sexes[[sx]]
  # 
  #   # if (dataset == "clinical") {
  #   #   # setting up background rects
  #   #   data_min <- min(comb_df$mean_parts - comb_df$std_error) # comb with error makes sure rects cover error bars too
  #   #   data_max <- max(comb_df$mean_parts + comb_df$std_error)
  #   #   # data_min <- data_mins[[anl]]
  #   #   # data_max <- data_maxs[[anl]]
  #   #   if (data_min >= sex_ranges[[sx]][2] || data_max <= sex_ranges[[sx]][1]) { # case of no normal_range
  #   #     ranges <- data.frame(ystart = data_min,
  #   #                          yend = data_max,
  #   #                          col = "outside_normal")
  #   #   } else if (data_min >= sex_ranges[[sx]][1] && data_max <= sex_ranges[[sx]][2]) { # case of normal_range only
  #   #     ranges <- data.frame(ystart = data_min,
  #   #                          yend = data_max,
  #   #                          col = "normal_range")
  #   #   } else if (data_min >= sex_ranges[[sx]][1]) { # case of normal_range + outside_normal upper
  #   #     ranges <- data.frame(ystart = c(data_min, sex_ranges[[sx]][2]),
  #   #                          yend = c(sex_ranges[[sx]][2], data_max),
  #   #                          col = c("normal_range", "outside_normal"))
  #   #   } else if (data_max <= sex_ranges[[sx]][2]) { # case of normal_range + outside_normal lower
  #   #     ranges <- data.frame(ystart = c(data_min, sex_ranges[[sx]][1]),
  #   #                          yend = c(sex_ranges[[sx]][1], data_max),
  #   #                          col = c("outside_normal", "normal_range"))
  #   #   } else { # case normal_range + outside_normal on both sides
  #   #     ranges <- data.frame(ystart = c(data_min, sex_ranges[[sx]][1], sex_ranges[[sx]][2]),
  #   #                          yend = c(sex_ranges[[sx]][1], sex_ranges[[sx]][2], data_max),
  #   #                          col = c("outside_normal", "normal_range", "outside_normal"))
  #   #   }
  #   # }
  # 
  #   pdf(file.path(graph_dir, paste(analy, dataset, sx, "Comb_Avg.pdf", sep = "_")), width = 9, height = 6)
  # 
  #   # normalizing with each fiber baseline
  #   # assumed sorted, so baseline is first timepoint of each analyte-fiber set, and fibers adjacent
  # 
  #   print(comb_df %>%
  #       ggplot() + #ggplot(position = pd) +
  #       # geom_rect(data = ranges, aes(ymin = ystart, ymax = yend, xmin = -Inf, xmax = Inf, fill = col), alpha = 0.4) +
  #       # scale_fill_manual(values = c(outside_normal = "#ffcccb", normal_range = "#c6ff95")) +
  #       geom_line(position = pd, aes(x = point, y = mean_parts, group = fiber, color = fiber)) +
  #       geom_errorbar(aes(x = point, ymin = mean_parts - std_error, ymax = mean_parts + std_error,
  #                         group = fiber, color = fiber),
  #                     width = .5, position = pd) +
  #       theme(
  #         axis.text.x = element_text(size=10, angle=45, hjust = 1)) +
  #       ylab(paste(analy_full, unit, sep = "")))
  #   dev.off()
  # }
}