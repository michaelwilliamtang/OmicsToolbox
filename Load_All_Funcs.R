library(tidyverse)
library(stringr)

# ran varvectors
fibers <- scan("Metadata/Fibers.tsv", character(), quote = '', sep = "\t")
ids <- scan("Metadata/Ids.tsv", character(), quote = '', sep = "\t")
timepoints <- scan("Metadata/Timepoints.tsv", character(), quote = '', sep = "\t")
datasets <- scan("Metadata/Datasets.tsv", character(), quote = '', sep = "\t")
small_datasets <- scan("Metadata/Small_Datasets.tsv", character(), quote = '', sep = "\t")

# loads all core functions (labeled _FUNC), in working dir
wd_files <- list.files()
wd_funcs <- wd_files[grep("FUNC[.]R$", wd_files)]

# known issue: if one of funcs throws error, will stop and skip sourcing others
for (func in wd_funcs) source(func)

# we also conveniently load in sel bile acids, sel chols, sel_bac, sel_all_bac
load("Data/Normalized/Normalized_Log_LCInulin_metabolomics_df.RData")
all_bile_acids <- grep('[C|c]holic', colnames(normalized_metabolomics_df), value=TRUE);
sel_bile_acids <- scan("Select_Bile_Acids.tsv", character(), quote = '', sep = "\t")
sel_chols <- c("Cholesterol..Total", "LDL..Calculated.", "LDL.HDL.Ratio", "Non.HDL.Chol..Calc")
sel_bac <- c("Bacteroides", "Oscillibacter", "Desulfovibrio", "Roseburia")
load("Data/Tidy_Full/Tidy_Full_Arabinoxylan_metaphlan_df.RData")
all_bac <- tidy_df$analyte %>% unique()
sel_all_bac <- all_bac[!grepl("[.]", all_bac) & !grepl("k__", all_bac)]