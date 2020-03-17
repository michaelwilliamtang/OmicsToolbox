library(tidyverse)
library(stringr)

# ran varvectors
fibers <- scan(file.path("Metadata", "Fibers.tsv"), character(), quote = '', sep = "\t", quiet = T)
ids <- scan(file.path("Metadata", "Ids.tsv"), character(), quote = '', sep = "\t", quiet = T)
timepoints <- scan(file.path("Metadata", "Timepoints.tsv"), character(), quote = '', sep = "\t", quiet = T)
datasets <- scan(file.path("Metadata", "Datasets.tsv"), character(), quote = '', sep = "\t", quiet = T)
small_datasets <- scan(file.path("Metadata", "Small_Datasets.tsv"), character(), quote = '', sep = "\t", quiet = T)

# loads all core functions (labeled _FUNC), in working dir
wd_files <- list.files()
wd_funcs <- wd_files[grep("FUNC[.]R$", wd_files)]

# known issue: if one of funcs throws error, will stop and skip sourcing others
for (func in wd_funcs) source(func)

# load convenient analyte lists
source_dir <- file.path("Metadata", "Analyte_Lists")
all_bile_acids <- scan(file.path("Metadata", "Analyte_Lists", "All_Bile_Acids.tsv"), character(), quote = '', sep = "\t", quiet = T)
sel_bile_acids <- scan(file.path("Metadata", "Analyte_Lists", "Select_Bile_Acids.tsv"), character(), quote = '', sep = "\t", quiet = T)
sel_chols <- scan(file.path("Metadata", "Analyte_Lists", "Select_Chols.tsv"), character(), quote = '', sep = "\t", quiet = T)
sel_bac <- scan(file.path("Metadata", "Analyte_Lists", "Select_Bac.tsv"), character(), quote = '', sep = "\t", quiet = T)
all_bac <- scan(file.path("Metadata", "Analyte_Lists", "All_Bac.tsv"), character(), quote = '', sep = "\t", quiet = T)
sel_all_bac <- scan(file.path("Metadata", "Analyte_Lists", "Select_All_Bac.tsv"), character(), quote = '', sep = "\t", quiet = T)