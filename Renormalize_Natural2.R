# renormalizes a tidied/gathered @tmp_df with cols id, analyte, point, val to create col renorm_val
# safer version, accounts for missing baselines by simply giving an error of sorts
# also makes the renormalization centered at participant mean instead of 0
# THIS VER: vectorized for speed, esp. for large datasets like genef

renormalize_natural <- function(tmp_df) {
  require(tidyverse)
  
  # sort, so baseline is first timepoint of each analyte-fiber set, and fibers adjacent
  tmp_df <- tmp_df %>% arrange(analyte, id, point)
  
  # getting means per analyte and attaching
  part_mean_df <- tmp_df %>% filter(point == "Baseline") %>%
    group_by(analyte) %>%
    dplyr::summarise(avg_bas = mean(val)) %>%
    ungroup()
  part_mean_vec <- part_mean_df$avg_bas
  names(part_mean_vec) <- part_mean_df$analyte
  tmp_df$bas <- part_mean_vec[tmp_df$analyte]
  
  # isolate ids with missing baselines to add after merge -> not needed, better implementation used in merge later
  # tmp_ids <- tmp_df %>% filter(point == "Baseline") %>%
  #   select(id) %>%
  #   unique()
  # miss_ids <- !(ids %in% tmp_ids)
  # miss_df <- tmp_df %>% filter(point)
  
  # getting means per analyte per id and attaching
  part_mean_df2 <- tmp_df %>% filter(point == "Baseline") %>%
    group_by(analyte, id) %>%
    dplyr::summarise(avg_bas2 = mean(val)) %>%
    ungroup()
  # nuance: for missing baselines the row for that id is omitted
  tmp_df2 <- merge(tmp_df, part_mean_df2, all.x = T) %>%
    arrange(analyte, id, point) # easier for visual comparison, optional
  print(table(tmp_df2$val == tmp_df$val)) # check merge
  
  # OLD RENORM, NON-VECTORIZED
  # subtract
  # bas <- 0
  # id <- ""
  # tmp_df <- tmp_df %>%
  #   group_by(id) %>%
  # for (i in 1:nrow(tmp_df)) {
  #   if (tmp_df[i,"point"] == "Baseline" && id != tmp_df[i, "id"]) {
  #     bas <- tmp_df[i, "val"]
  #     id <- tmp_df[i, "id"]
  #     # print(id)
  #   }
  #   # print(bas)
  #   if (tmp_df[i, "id"] == id) {
  #     tmp_df[i, "renorm_val"] <- tmp_df[i, "val"] - bas
  #   } else {
  #     print(paste("WARNING MISSING BASELINE", paste(tmp_df[i, ], collapse = " ")))
  #   }
  # }
  
  # norm and recentering to part means
  tmp_df2$renorm_val <- tmp_df2$val - tmp_df2$avg_bas2+ tmp_df2$bas
  tmp_df2 <- tmp_df2 %>% select(-c(bas, avg_bas2)) # no need after norm done
  
  return(tmp_df2)
}