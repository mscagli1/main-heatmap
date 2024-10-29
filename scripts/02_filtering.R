# Purpose, Input and Outputs ----------------------------------------------
# The purpose of this script is to filter the main data frame for DEGs.
# Input: main data frame from script01
# Output: R object: mysample_groups - list of deprivation or drug condition vectors
#         R objectL control_CPMs - list of control conditions for retaining mean count data
#         R object: RNAdata_deps.DE.rds list of dataframes containing differential genes by condition and FC cutoff
#         .csv files of DEGS for each condition

# Load Packages and set directory -----------------------------------------
here::i_am("scripts/02_filtering.R")
library(here) #Creates directory-specific paths
library(tidyverse)

# Load Data ---------------------------------------------------------------
main_data <- readRDS(here("output", "objects", "s01_main_data.rds")) #Load in renamed data

# Create strings and testing columns --------------------------------------
# Create strings representing variable names for later iteration and filtering
sel_deps      <- readRDS(here("output", "objects", "s01_sel_deps.rds"))
type          <- readRDS(here("output", "objects", "s01_type.rds"))

# test_columns is a list, containing lists of column name vectors, grouped by sample group (sel_deps/drugs/all)
# Each group contains paired vectors of condition names, FC column names, and pvalue column names for iterating on!
# The purpose of this is to select genes that are significant and differential for the same condition!
test_columns         <- list()

create_types <- function(conditions) {
  new_vec  <- c(paste0(conditions, ".TotalRNA"),
                paste0(conditions, ".PolyRNA"))
  return(new_vec)
}

test_columns$sample_group <- sel_deps |> map(create_types) |> unlist()
test_columns$FC_cols      <- test_columns$sample_group |> paste0(".Log2FC")
test_columns$adjP_cols    <- test_columns$sample_group |> paste0(".AdjPval")

saveRDS(test_columns, here("output", "objects", "s02_test_columns.rds"))

# Filtering ---------------------------------------------------------------
#Select RNA data columns and filter out NA or duplicated gene name data
RNAdata_deps <- main_data |>
  select(contains(c("Name", "biotype", "RNA"))) |>
  select(!matches(c("rawCounts", "normCounts", "CPM"), ignore.case = FALSE)) |>
  filter(!is.na(Name)) |>
  group_by(.data$Name) |>
  filter(n() == 1) |>
  ungroup()

#Select only data columns in data group (sel_deps, drugs, etc) and filter for detected genes across every sample_group
RNAdata_deps <- RNAdata_deps |>
    select(contains(c("Name", sel_deps, "Control"))) |>
    filter(if_all(contains("Log2FC"), ~ !is.na(.x)))

saveRDS(RNAdata_deps, here("output", "objects", "s02_RNAdata_deps.rds"))
write_csv(RNAdata_deps, here("data", "processed", "s02_RNAdata_deps.csv"))

# FoldChange and pValue cutoffs for filtering
my_FC_cutoff <- 1
my_padj_cutoff <- 0.05

#Differential gene testing
tempDEGenes <- list()
  # For each data group, create a list of dataframes, with each data frame consisting of all DE/sig genes from that condition
  # Bind all dataframes together, keeping only distinct rows to avoid duplicate data.
  # This allows us to make discrete comparisons of one FC column AND one pval column for the same condition

test_DEGs <- function (data, sample_group, FC_cols, adjP_cols, FC_cutoff, padj_cutoff) {
      df <- data |>
    filter((.data[[FC_cols]] > FC_cutoff | .data[[FC_cols]] < -FC_cutoff) & .data[[adjP_cols]] <= padj_cutoff)
  return(df)
}

tempDEGenes <- pmap(data = RNAdata_deps,
                    .l = test_columns,
                    .f = test_DEGs,
                    FC_cutoff = my_FC_cutoff,
                    padj_cutoff = my_padj_cutoff
                    )
tempDEGenes <- set_names(tempDEGenes, test_columns$sample_group)


RNAdata_deps_DE <- list()
RNAdata_deps_DE <- bind_rows(tempDEGenes) |> distinct()  # Keep only distinct rows and overwrite

saveRDS(RNAdata_deps_DE, here("output", "objects", "s02_RNAdata_deps_DE.rds"))
write_csv(RNAdata_deps_DE, here("data", "processed", "s02_RNAdata_deps_DE.csv"))

message("You did it!")
