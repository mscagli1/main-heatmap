# Purpose, Input and Output -----------------------------------------------
  # Prepare data for Log2FC heatmapping and pathway-specific heatmapping

# Load packages -----------------------------------------------------------
here::i_am("scripts/05_ko-heatmaps.R")
library(here) #Creates directory-specific paths
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(gprofiler2)
library(qpdf)
library(quarto)

# Load data and objects ---------------------------------------------------
ko_only_DE      <- readRDS(here("output", "objects", "s04_ko_only_DE.rds"))
ko_halo_DE      <- readRDS(here("output", "objects", "s04_ko_halo_DE.rds"))
test_columns_ko_halo <- readRDS(here("output", "objects", "s04_test_columns_ko_halo.rds"))
test_columns_ko_only <- readRDS(here("output", "objects", "s04_test_columns_ko_only.rds"))

# pathway_genes   <- readRDS(here("data", "input", "genes.in.gene.sets.of.interest.rds"))
# pathway_annos   <- read_csv(here("data", "input", "my_selected_gene_sets_mapped.csv"), col_select = (c(database, my_name, id)))



# Split data --------------------------------------------------------------

#Split data into count and fc data
control_cpms <- c("DMSO.Rosa26.meanNormCounts", "Halo.Rosa26.meanNormCounts")
ko_only_DE_FC     <- ko_only_DE |> select(contains(c("Name", "Log2FC", control_cpms)))
ko_only_DE_counts <- ko_only_DE |> select(contains(c("Name", "meanNormCounts")))

ko_halo_DE_FC     <- ko_halo_DE |> select(contains(c("Name", "Log2FC", control_cpms)))
ko_halo_DE_counts <- ko_halo_DE |> select(contains(c("Name", "meanNormCounts")))

# Prep data ------------------------------------------------------------

# Preps data for heatmapping by selecting columns, filtering FC,
# ordering columns, splitting FC and CPM data into separate objects,
# and returning FC and CPM data in a list.
prep_hm <- function(data) {

  filt_data <- data |> filter(if_all(contains("meanNormCounts"), ~ .x > 0.5))

  FC_data <- filt_data |>
    select(contains(c("Name", "ATF4-vs-Rosa26", "CEBPG-vs-Rosa26"))) |>
    relocate(contains(c("DMSO", "Halo"))) |>
    column_to_rownames("Name")

  CPM_data <- filt_data  |>
    select(contains(c("Name", "meanNormCounts"))) |>
    column_to_rownames("Name")

  Halo_data <- filt_data  |>
    select(contains(c("Name", "Halo-vs-DMSO"))) |>
    column_to_rownames("Name")

  output <- list(FCs = FC_data, control_CPMs = CPM_data, Halo = Halo_data)
  return(output)
}

ko_only_DE_FC_hm <- ko_only_DE_FC |> prep_hm()
ko_halo_DE_FC_hm <- ko_halo_DE_FC |> prep_hm()

prep_hm_halo <- function(data) {

  filt_data <- data |> filter(if_all(contains("meanNormCounts"), ~ .x > 0.5))

  FC_data <- filt_data |>
    select(contains(c("Name", "Log2FC"))) |>
    relocate(contains(c("DMSO", "Halo"))) |>
    column_to_rownames("Name")

  CPM_data <- filt_data  |>
    select(contains(c("Name", "meanNormCounts"))) |>
    column_to_rownames("Name")

  output <- list(FCs = FC_data, control_CPMs = CPM_data)
  return(output)
}

ko_halo2_DE_FC_hm <- ko_halo_DE_FC |> prep_hm_halo()

# Clustering --------------------------------------------------------------

# This function takes data and a given number of clusters (k). It returns a set of cluster assignments
# and ranks them to make cluster order and name more reproducible as k changes.
cluster_FC_data <- function(data, k) {
  set.seed(201) #Set seed for more reproducible results
  myclust <- kmeans(data$FCs, k, iter.max = 500) #Set data, number of clusters, number of iterations
  myclust$k <- k

  # Rename clusters in order of highest sum of FCs to lowest
  order <- names(sort(rowSums(myclust$centers), decreasing = TRUE))
  names(order) <- as.character(seq_along(order))

  # Reassigns cluster names in decreasing order of sum
  names <- names(myclust$cluster)
  old_cluster <- myclust$cluster
  new_cluster <- fct_recode(as.character(old_cluster), !!!order)
  names(new_cluster) <- names

  # Saves new cluster names as a factor in the clust object
  ord_clust <- factor(new_cluster, levels = seq_along(order))
  myclust$ord_clust <- ord_clust

  # Save gene names for each cluster
  myclust$genes_by_ord_clust <- list()
  for (i in 1:k) {
    cluster_name <- paste0("Cluster", sep = "_", i)
    myclust$genes_by_ord_clust[[cluster_name]] <- names(ord_clust[ord_clust == i])
  }

  return(myclust)
}

# Cluster my data
ko_only_DE_FC_hm$clust <- ko_only_DE_FC_hm |> cluster_FC_data(4)
ko_halo_DE_FC_hm$clust <- ko_halo_DE_FC_hm |> cluster_FC_data(8)

ko_halo2_DE_FC_hm$clust <- ko_halo2_DE_FC_hm |> cluster_FC_data(8)



# Make global Heatmap with halo data as annotation -----------------------------------------------------

createKO_Heatmap <- function(data, filename) {

    fc_matrix <- as.matrix(data$FCs)
    cpm_matrix <- as.matrix(data$control_CPMs)
    halo_matrix <- as.matrix(data$Halo)
    cond_labels <- c("ATF4 KO vs. Rosa26 (in DMSO)",
                    "CEBPG KO vs. Rosa26 (in DMSO)",
                    "ATF4 KO vs. Rosa26 (in Halo)",
                    "CEBPG KO vs. Rosa26 (in Halo)"
    )
    cpm_labels <- c("DMSO (Rosa26)", "Halo (Rosa26)")
    colnames(cpm_matrix) <- cpm_labels

    halo_labels <- "Log2FC - Halo vs DMSO (in Rosa26)"
    colnames(halo_matrix) <- halo_labels

    k <- data$clust$k

    cond_colors <- c("red2", "dodgerblue2", "red","dodgerblue2")

    #Main data color mapping
    colors_main    <- colorRamp2(c(-2, 0, 2), c("dodgerblue2", "white", "red1"))
    #Right anno color mapping
    colors_CPM <- colorRamp2(c(0, 4, 8, 15), c("white", "gold1", "orange", "red1"))

    colors_cluster <- brewer.pal(n = k, name = "Dark2")
    names(colors_cluster) <- as.character(1:k)

    colors_halo    <- colorRamp2(c(-2, 0, 2), c("dodgerblue2", "white", "red1"))

    #Build final row annotations
    left_annos <- rowAnnotation(
      Cluster = data$clust$ord_clust,
      Expression         = cpm_matrix,
      Halo = halo_matrix,
      #Color assignments
      col         = list(Expression = colors_CPM,
                         Cluster    = colors_cluster,
                         Halo = colors_halo
      ),
      border = c(Expression = TRUE),
      annotation_name_offset = unit(2, "mm"),
      annotation_name_gp = gpar(fontsize = 10),
      annotation_width = c(1.5, 5, 10), width = unit(20, "mm"),
      simple_anno_size_adjust = TRUE,
      show_annotation_name = c(Cluster = FALSE)
    )


    selected_genes <- list()
    selected_genes$C1 <- c("Slc7a11", "Slc3a2")
    selected_genes$C2 <- c("Gadd45g")
    selected_genes$C3 <- c("Asns", "Bcat1", "Shmt2", "Atf4", "Cebpg", "Cars1")
    selected_genes$C4 <- c("Fosl1", "Junb", "Mthfd1l", "Lag3", "Rars1", "Bub1b")
    selected_genes$C5 <- c("Irf4", "Socs3", "Slc2a1")
    selected_genes$C6 <- c("Rps25", "Rpl35", "Rpl32", "Mrps36", "Eif1b", "Ndufs6","Ndufa4")
    selected_genes$C7 <- c("Slc16a3", "Bckdk", "Fabp5", "Il2ra")
    selected_genes$C8 <- c("Hmgcs1")


    my_favorite_genes <- unlist(selected_genes, use.names = FALSE)

    marked_indices <- which((rownames(fc_matrix) %in% my_favorite_genes) == "TRUE")
    marked_labels <- rownames(fc_matrix)[marked_indices]

    right_annos <- rowAnnotation(
      Genes = anno_mark(at = marked_indices,
                        labels    = marked_labels,
                        labels_gp = gpar(fontsize = 10)),
      annotation_name_offset = unit(2, "mm"),
      width = unit(10, "mm"),
      show_annotation_name = c(Genes = FALSE)
    )

    splits <- data$clust$ord_clust

    ht_opt$TITLE_PADDING = unit(4, "mm")
    ht_opt$DIMNAME_PADDING = unit(2, "mm")
    ht_opt$ROW_ANNO_PADDING = unit(2, "mm")
    heatmap <- Heatmap(matrix            = fc_matrix,
                       col                = colors_main,
                       left_annotation    = left_annos,
                       right_annotation = right_annos,
                       name = str_wrap("log2FC vs Control", width = 10),

                       heatmap_width      = unit(9, "cm"),
                       heatmap_height     = unit(16, "cm"),

                       column_order       = colnames(data$FCs),
                       row_split = splits,
                       row_title_rot = 45,
                       column_labels = cond_labels,
                       show_row_dend      = FALSE,
                       show_row_names     = FALSE,
                       row_gap            = unit(1, "mm"),
                       column_gap         = unit (2, "mm"),
                       column_names_gp    = gpar(col = cond_colors, fontface = "bold"),
                       row_names_gp       = gpar(fontsize = 8),
                       border             = TRUE,

                       cluster_columns    = FALSE,
                       cluster_row_slices = FALSE,
                       use_raster         = TRUE
    )

    plot <- draw(heatmap, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")

    saveRDS(plot, here("output", "objects", paste0(filename, ".rds")))

    pdf(here("output", "plots", paste0(filename,".pdf")), width = 5, height = 9)
    plot(plot)
    dev.off()

    ht_opt(RESET =  TRUE)
    return(plot)
  }


ko_only_heatmap <- createKO_Heatmap(ko_only_DE_FC_hm, filename = "s05_ko_only_heatmap")
ko_halo_heatmap <- createKO_Heatmap(ko_halo_DE_FC_hm, filename = "s05_ko_halo_heatmap")

createKO_Heatmap2 <- function(data, filename) {

  fc_matrix <- as.matrix(data$FCs)
  cpm_matrix <- as.matrix(data$control_CPMs)
  cond_labels <- c("Halo vs DMSO (Rosa26)",
                    "ATF4 KO vs. Rosa26 (in DMSO)",
                   "CEBPG KO vs. Rosa26 (in DMSO)",
                   "ATF4 KO vs. Rosa26 (in Halo)",
                   "CEBPG KO vs. Rosa26 (in Halo)"
  )
  cpm_labels <- c("DMSO (Rosa26)", "Halo (Rosa26)")
  colnames(cpm_matrix) <- cpm_labels

  k <- data$clust$k

  cond_colors <- c("springgreen3", "red2", "dodgerblue2", "red","dodgerblue2")

  #Main data color mapping
  colors_main    <- colorRamp2(c(-2, 0, 2), c("dodgerblue2", "white", "red1"))
  #Right anno color mapping
  colors_CPM <- colorRamp2(c(0, 4, 8, 15), c("white", "gold1", "orange", "red1"))

  colors_cluster <- brewer.pal(n = k, name = "Dark2")
  names(colors_cluster) <- as.character(1:k)

  #Build final row annotations
  left_annos <- rowAnnotation(
    Cluster = data$clust$ord_clust,
    Expression         = cpm_matrix,
    #Color assignments
    col         = list(Expression = colors_CPM,
                       Cluster    = colors_cluster
    ),
    border = c(Expression = TRUE),
    annotation_name_offset = unit(2, "mm"),
    annotation_name_gp = gpar(fontsize = 10),
    annotation_width = c(1.5, 5), width = unit(10, "mm"),
    simple_anno_size_adjust = TRUE,
    show_annotation_name = c(Cluster = FALSE)
  )


  selected_genes <- list()
  selected_genes$C1 <- c("Slc7a11", "Slc3a2")
  selected_genes$C2 <- c("Gadd45g")
  selected_genes$C3 <- c("Asns", "Bcat1", "Shmt2", "Atf4", "Cebpg", "Cars1")
  selected_genes$C4 <- c("Fosl1", "Junb", "Mthfd1l", "Lag3", "Rars1", "Bub1b")
  selected_genes$C5 <- c("Irf4", "Socs3", "Slc2a1")
  selected_genes$C6 <- c("Rps25", "Rpl35", "Rpl32", "Mrps36", "Eif1b", "Ndufs6","Ndufa4")
  selected_genes$C7 <- c("Slc16a3", "Bckdk", "Fabp5", "Il2ra")
  selected_genes$C8 <- c("Hmgcs1")


  my_favorite_genes <- unlist(selected_genes, use.names = FALSE)

  marked_indices <- which((rownames(fc_matrix) %in% my_favorite_genes) == "TRUE")
  marked_labels <- rownames(fc_matrix)[marked_indices]

  right_annos <- rowAnnotation(
    Genes = anno_mark(at = marked_indices,
                      labels    = marked_labels,
                      labels_gp = gpar(fontsize = 10)),
    annotation_name_offset = unit(2, "mm"),
    width = unit(10, "mm"),
    show_annotation_name = c(Genes = FALSE)
  )

  splits <- data$clust$ord_clust

  ht_opt$TITLE_PADDING = unit(4, "mm")
  ht_opt$DIMNAME_PADDING = unit(2, "mm")
  ht_opt$ROW_ANNO_PADDING = unit(2, "mm")
  heatmap <- Heatmap(matrix            = fc_matrix,
                     col                = colors_main,
                     left_annotation    = left_annos,
                     right_annotation = right_annos,
                     name = str_wrap("log2FC vs Control", width = 10),

                     heatmap_width      = unit(9, "cm"),
                     heatmap_height     = unit(16, "cm"),

                     column_order       = colnames(data$FCs),
                     row_split = splits,
                     row_title_rot = 45,
                     column_labels = cond_labels,
                     show_row_dend      = FALSE,
                     show_row_names     = FALSE,
                     row_gap            = unit(1, "mm"),
                     column_gap         = unit (2, "mm"),
                     column_names_gp    = gpar(col = cond_colors, fontface = "bold"),
                     row_names_gp       = gpar(fontsize = 8),
                     border             = TRUE,

                     cluster_columns    = FALSE,
                     cluster_row_slices = FALSE,
                     use_raster         = TRUE
  )

  plot <- draw(heatmap, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")

  saveRDS(plot, here("output", "objects", paste0(filename, ".rds")))

  pdf(here("output", "plots", paste0(filename,".pdf")), width = 5, height = 9)
  plot(plot)
  dev.off()

  ht_opt(RESET =  TRUE)
  return(plot)
}


ko_halo2_heatmap <- createKO_Heatmap2(ko_halo2_DE_FC_hm, filename = "s05_ko_halo2_heatmap")


# Cluster gProfiler Plots -----------------------------------------------------------

plot_gostplot <- function(genes) {

  gost_result <- gost(query = genes,
                      organism = "mmusculus", ordered_query = FALSE,
                      multi_query = FALSE, significant = TRUE, exclude_iea = TRUE,
                      measure_underrepresentation = FALSE, evcodes = FALSE,
                      user_threshold = 0.05, correction_method = "g_SCS",
                      domain_scope = "annotated", custom_bg = NULL,
                      sources = NULL, as_short_link = FALSE, highlight = FALSE)

  gost_plot <- gostplot(gost_result, capped = TRUE, interactive = TRUE)
  return(gost_plot)
}


plot_gostplots_by_clusters <- function(clust_gene_list, filename) {

  my_plots <- vector("list", length(clust_gene_list))
  names(my_plots) <- names(clust_gene_list)

  for (name in names(clust_gene_list)) {
    my_plots[[name]] <- plot_gostplot(clust_gene_list[name])

  }
    saveRDS(my_plots, here("output", "objects", paste0(filename, ".rds")))
    saveRDS(clust_gene_list, here("output", "objects", paste0(filename, "_clust_genes.rds")))

    return(my_plots)
}

gostplots_ko_halo <- plot_gostplots_by_clusters(ko_halo_DE_FC_hm$clust$genes_by_ord_clust, "s05_ko_halo_gostplots")
gostplots_ko_only <- plot_gostplots_by_clusters(ko_only_DE_FC_hm$clust$genes_by_ord_clust, "s05_ko_only_gostplots")
gostplots_ko_halo2 <- plot_gostplots_by_clusters(ko_halo2_DE_FC_hm$clust$genes_by_ord_clust, "s05_ko_halo2_gostplots")





# Render Report -----------------------------------------------------------


quarto_render(
  input = here("output", "reports", "2024-10-22_ko-heatmap-clusters.qmd"),
)
