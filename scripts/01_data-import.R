# Purpose, inputs, and outputs -------------------------------------------------------------------

# The purpose of this script is to load in Montana's main data frame for the RNA-seq analysis (deprivation polysome data, drug polysome data, and KO RNA-seq),
# alter column names for easy use with tidy-select conventions, and create derivative datasets for downstream analysis.

# Input: Montana's "main" data frame
# Outputs:
#   column mapping to trace names of old vs new columns
#   renamed main data file to split into new data frames downstream

# Load packages and set "here" directory ----------------------------------
here::i_am("scripts/01_data-import.R")
library(here) #Creates directory-specific paths
library(tidyverse)

# Load data ---------------------------------------------------------------
main_data_path <- here("data", "input", "Deprivation.and.DrugTreatment.and.Knockout.LogFC.and.Anota2seqStatuses.and.Counts.csv")
main_data      <- read_csv(main_data_path)

# Edit column name strings --------------------------------------------------------
all_colnames_old <- colnames(main_data)

colnames(main_data)[-1:-3] <-  colnames(main_data[-1:-3]) |>
  str_replace_all(c(" " = ".",
                    "_" = ".",
                    "-" = "."
                    )
                  )

colnames(main_data)[-1:-3] <- colnames(main_data)[-1:-3] |>
  str_replace_all(c("\\.adjusted\\.pval"     = ".AdjPval",
                    "\\.pval"                = ".AdjPval",
                    "\\.Pval"                = ".AdjPval",
                    "\\.AdjPvalue"           = ".AdjPval",

                    "logFC"                  = "Log2FC",
                    "LogFC"                  = "Log2FC",
                    "logCPM"                 = "LogCPM",
                    "Anota2seq\\.Status"     = "Anota2seqStatus",
                    "\\.EGi1"                = ".FourEGi1",
                    "4EGi1"                  = "FourEGi1",
                    "Halofuginone"           = "Halo",
                    "^treated"               = "Halo",
                    "Ctrl"                   = "DMSO",
                    "untreated"              = "DMSO",

                    "Protein\\.NoGlucose"    = "NoGlucose\\.Protein",
                    "Protein\\.NoGlutamine"  = "NoGlutamine\\.Protein",
                    "Protein\\.NoBCAA"       = "NoBCAA\\.Protein",
                    "Protein\\.NoArginine"   = "NoArginine\\.Protein",
                    "Protein\\.NoMethionine" = "NoMethionine\\.Protein",
                    "Protein\\.Torin"        = "Torin\\.Protein",
                    "Protein\\.Halo"         = "Halo\\.Protein",
                    "Protein\\.FourEGi1"     = "FourEGi1\\.Protein",
                    "Protein\\.S6Ki"         = "S6Ki\\.Protein",


                    "IN"                     = "TotalRNA",
                    "PO"                     = "PolyRNA",
                    "Input"                  = "TotalRNA",
                    "PooledPolysome"         = "PolyRNA",
                    "raw\\.counts"           = "rawCounts",
                    "normalized\\.counts"    = "normCounts",
                    "rosa26"                 = "Rosa26",
                    "atf4"                   = "ATF4",
                    "cebpg"                  = "CEBPG",
                    "Rosa\\."                = "Rosa26\\.",
                    "^[0-9][0-9]\\."         = "",

                    "Halo.vs.DMSO" = "Halo-vs-DMSO",
                    "ATF4.vs.Rosa26" = "ATF4-vs-Rosa26",
                    "CEBPG.vs.Rosa26" = "CEBPG-vs-Rosa26"
                    )
                  )

#Repair Jun data
jun_index <- which(main_data$ensembl_gene_id_version == "ENSMUSG00000052684.5")
main_data$Name[jun_index] <- "Jun"
main_data$gene_biotype[jun_index] <- "protein_coding"

# Repair Rep info ----------------------------------------------------------------
#Rename normalized count columns with missing rep data
missingNormCounts_names <-colnames(main_data[141:176]) #These columns have missing rep data.

rep_info <- c(rep("Rep1", 6), rep("Rep2", 6), rep("Rep3", 6),
              rep("Rep1", 6), rep("Rep2", 6), rep("Rep3", 6)
              )

missingNormCounts_newnames <- vector()
for (i in 1:length(missingNormCounts_names)){
  info <- missingNormCounts_names |> str_split("\\.", simplify = FALSE)
  type <- info[[i]][[2]]
  cond <- info[[i]][[3]]
  measure <- info[[i]][[4]]
  rep <- rep_info[i]
  missingNormCounts_newnames[i] <- paste(type, cond, rep, measure, sep = ".")
}

colnames(main_data)[141:176] <- missingNormCounts_newnames #Apply new names
all_colnames_new <- colnames(main_data) #Save all names


# Swap order of count data names --------------------------------------
all_deps          <- c("Control", "NoGlucose", "NoGlutamine", "NoBCAA", "NoArginine", "NoMethionine")
all_drugs         <- c("DMSO", "Torin", "Halo", "FourEGi1")
allconditions <- c(all_deps, all_drugs)

sel_deps          <- c("NoGlucose", "NoGlutamine", "NoBCAA", "NoArginine", "NoMethionine")
sel_drugs         <- c("Torin", "Halo", "FourEGi1")

type          <- c("TotalRNA", "PolyRNA")

#Save objects
saveRDS(all_deps, here("output", "objects", "s01_deps.rds"))
saveRDS(all_drugs, here("output", "objects", "s01_drugs.rds"))
saveRDS(allconditions, here("output", "objects", "s01_allconditions.rds"))
saveRDS(sel_deps, here("output", "objects", "s01_sel_deps.rds"))
saveRDS(sel_drugs, here("output", "objects", "s01_sel_drugs.rds"))
saveRDS(type, here("output", "objects", "s01_type.rds"))


switched_cols_old <- colnames(main_data) |>
  str_subset("normCounts$|rawCounts$") |>
  str_subset(paste(type, collapse = "|")) |>
  str_subset(paste(allconditions, collapse = "|"))

switched_cols_new <-
  paste(str_split_i(switched_cols_old, "\\.", 2),
        str_split_i(switched_cols_old, "\\.", 1),
        str_split_i(switched_cols_old, "\\.", 3),
        str_split_i(switched_cols_old, "\\.", 4),
        sep = ".")

cols_to_switch <- which(colnames(main_data) %in% switched_cols_old)
colnames(main_data)[cols_to_switch] <- switched_cols_new

# Create mapping and save ----------------------------------------------------------
#Create column mapping in case Montana has the original names in the dataframe
column_mapping <- data.frame(old = all_colnames_old, new = all_colnames_new)
write.csv(column_mapping, here("data", "processed", "s01_main_data_column_mapping.csv"), row.names = FALSE)


# Create mean of normalized counts ----------------------------------------
norm_count_cols <- colnames(main_data) |> str_subset("normCounts$")
norm_count_groups <- paste(str_split_i(norm_count_cols, "\\.", 1), str_split_i(norm_count_cols, "\\.", 2), sep = ".") |>
  unique()

for (name in norm_count_groups){
  main_data <- main_data |>
    mutate("{name}.meanNormCounts" := rowMeans(pick(contains("normCounts") & matches(name))))
  }


#Save data
write.csv(main_data, here("data", "processed", "s01_main_data.csv"), row.names = FALSE)
saveRDS(main_data, here("output", "objects", "s01_main_data.rds"))

message("Data import - Done!")

