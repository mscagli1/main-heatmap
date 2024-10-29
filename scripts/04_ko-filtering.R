# Purpose, Input and Outputs ----------------------------------------------
# The purpose of this script is to filter the main data frame for DEGs.
# Input: main data frame from script01
# Output: R object: mysample_groups - list of deprivation or drug condition vectors
#         R objectL control_CPMs - list of control conditions for retaining mean count data
#         R object: ko_data.DE.rds list of dataframes containing differential genes by condition and FC cutoff
#         .csv files of DEGS for each condition

# Load Packages and set directory -----------------------------------------
here::i_am("scripts/04_ko-filtering.R")
library(here) #Creates directory-specific paths
library(tidyverse)

# Script Number -----------------------------------------------------------
scr_no <- "s04"

# Load Data ---------------------------------------------------------------
main_data <- readRDS(here("output", "objects", "s01_main_data.rds")) #Load in renamed data

ko_data <- main_data |>
  select(contains(c("Name", "Rosa26", "ATF4", "CEBPG"))) |>
  select(contains(c("Name", "Log2FC", "AdjPval", "meanNormCounts"))) |>
  select(!contains(c("Halo-vs-DMSO.ATF4-vs-Rosa26",
                     "Halo-vs-DMSO.CEBPG-vs-Rosa26",
                     "Halo-vs-DMSO.ATF4",
                     "Halo-vs-DMSO.CEBPG")))


#Select RNA data columns and filter out NA or duplicated gene name data
ko_data <- ko_data |>
  filter(!is.na(Name)) |>
  group_by(.data$Name) |>
  filter(n() == 1) |>
  ungroup() |>
  filter(if_all(contains("Log2FC"), ~ !is.na(.x)))

ko_data_counts  <- ko_data |> select(contains(c("Name", "meanNormCounts")))
ko_data_fc_p <- ko_data |> select(contains(c("Name", "Log2FC", "AdjPval")))

saveRDS(ko_data, here("output", "objects", paste0(scr_no, sep = "_", "ko_data.rds")))
write_csv(ko_data, here("data", "processed", paste0(scr_no, sep = "_", "ko_data.csv")))

# test_columns_ko_halo is a list, containing lists of column name vectors, grouped by sample group (sel_deps/drugs/all)
# Each group contains paired vectors of condition names, FC column names, and pvalue column names for iterating on!
# The purpose of this is to select genes that are significant and differential for the same condition!

test_columns_ko_halo         <- list()
test_columns_ko_halo$sample_group <- str_remove(names(ko_data_fc_p)[2:6], ".Log2FC")
test_columns_ko_halo$FC_cols      <- test_columns_ko_halo$sample_group |> paste0(".Log2FC")
test_columns_ko_halo$adjP_cols    <- test_columns_ko_halo$sample_group |> paste0(".AdjPval")

test_columns_ko_only         <- list()
test_columns_ko_only$sample_group <- str_remove(names(ko_data_fc_p)[3:6], ".Log2FC")
test_columns_ko_only$FC_cols      <- test_columns_ko_only$sample_group |> paste0(".Log2FC")
test_columns_ko_only$adjP_cols    <- test_columns_ko_only$sample_group |> paste0(".AdjPval")


check_test_cols <- function(data, cols) {
  test_vec <- cols %in% names(data)
  if (all(test_vec) == TRUE) {
    message("Test columns found!")
    } else {
      stop("Test columns not found - check your column names!")
    }
  }

check_test_cols(data = ko_data_fc_p, cols = test_columns_ko_halo$FC_cols)
check_test_cols(data = ko_data_fc_p, cols = test_columns_ko_halo$adjP_cols)
check_test_cols(data = ko_data_fc_p, cols = test_columns_ko_only$FC_cols)
check_test_cols(data = ko_data_fc_p, cols = test_columns_ko_only$adjP_cols)

saveRDS(test_columns_ko_halo, here("output", "objects", paste0(scr_no, sep = "_", "test_columns_ko_halo.rds")))
saveRDS(test_columns_ko_only, here("output", "objects", paste0(scr_no, sep = "_", "test_columns_ko_only.rds")))

# Differential gene testing ---------------------------------------------------------------

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



tempDEGenes <- list()
tempDEGenes <- pmap(data = ko_data,
                    .l = test_columns_ko_halo,
                    .f = test_DEGs,
                    FC_cutoff = my_FC_cutoff,
                    padj_cutoff = my_padj_cutoff
                    )
tempDEGenes <- set_names(tempDEGenes, test_columns_ko_halo$sample_group)

ko_halo_DE <- list()
ko_halo_DE <- bind_rows(tempDEGenes) |> distinct()  # Keep only distinct rows and overwrite



tempDEGenes2 <- list()
tempDEGenes2 <- pmap(data = ko_data,
                    .l = test_columns_ko_only,
                    .f = test_DEGs,
                    FC_cutoff = my_FC_cutoff,
                    padj_cutoff = my_padj_cutoff
)
tempDEGenes2 <- set_names(tempDEGenes2, test_columns_ko_only$sample_group)

ko_only_DE <- list()
ko_only_DE <- bind_rows(tempDEGenes2) |> distinct()  # Keep only distinct rows and overwrite



saveRDS(ko_halo_DE, here("output", "objects", paste0(scr_no, sep = "_", "ko_halo_DE.rds")))
write_csv(ko_halo_DE, here("data", "processed", paste0(scr_no, sep = "_", "ko_halo_DE.csv")))

saveRDS(ko_only_DE, here("output", "objects", paste0(scr_no, sep = "_", "ko_only_DE.rds")))
write_csv(ko_only_DE, here("data", "processed", paste0(scr_no, sep = "_", "ko_only_DE.csv")))

message("You did it!")
