# Purpose, Input and Output -----------------------------------------------
# Create a violin plot to show gene expression by biotype

# Load packages -----------------------------------------------------------
here::i_am("scripts/06_biotype-plot.R")
library(here) #Creates directory-specific paths
library(tidyverse)
library(ggplot2)
library(introdataviz)


# Import Data -------------------------------------------------------------
main_data <- readRDS(here("output", "objects", "s01_main_data.rds"))

mean_counts <- main_data |>
  select(contains(c("Name", "biotype", "meanNormCounts"))) |>
  select(contains(c("Name", "biotype", "Control"))) |>
  pivot_longer(cols = contains("Control"), names_to = "sample", values_to = "meanNormCounts") |>
  drop_na(meanNormCounts, Name, gene_biotype)


mean_counts$sample <- mean_counts$sample |>
  case_match("Control.TotalRNA.meanNormCounts" ~ "Total RNA",
             "Control.PolyRNA.meanNormCounts" ~ "Poly RNA")

sums <- mean_counts |>
  group_by(gene_biotype) |>
  summarize(count = n())

keep <- c("TEC",
          "lncRNA",
          "processed_pseudogene",
          "protein_coding",
          "snRNA",
          "snoRNA",
          "unprocessed_pseudogene")

mean_counts <- mean_counts |>
  filter(gene_biotype %in% keep)

mean_counts$gene_biotype <- mean_counts$gene_biotype |>
  str_replace(pattern = "_", replacement = " ")

biotype_plot <- mean_counts |>
  ggplot(aes(x = gene_biotype, y = meanNormCounts, fill = factor(sample, levels = c("Total RNA", "Poly RNA")))) +
  introdataviz::geom_split_violin(alpha = .4, trim = FALSE) +
  geom_boxplot(width = .2, alpha = .6, outliers = FALSE, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F,
               position = position_dodge(.175)) +
  scale_x_discrete(name = "RNA type") +
  scale_y_continuous(name = "logCPM") +
  scale_fill_brewer(palette = "Dark2", name = "RNA type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(here("output", "plots", "s06_biotype_plot.pdf"), height = 5, width = 7)
