# Purpose, Input and Output -----------------------------------------------
  # Prepare data for Log2FC heatmapping and pathway-specific heatmapping

# Load packages -----------------------------------------------------------
here::i_am("scripts/03_main-heatmap.R")
library(here) #Creates directory-specific paths
library(tidyverse)
library(glue)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(gprofiler2)
library(qpdf)

# Load data and objects ---------------------------------------------------
RNAdata_deps      <- readRDS(here("output", "objects", "s02_RNAdata_deps.rds"))
RNAdata_deps_DE   <- readRDS(here("output", "objects", "s02_RNAdata_deps_DE.rds"))
test_columns      <- readRDS(here("output", "objects", "s02_test_columns.rds"))
sel_deps          <- readRDS(here("output", "objects", "s01_sel_deps.rds"))
pathway_genes     <- readRDS(here("data", "input", "genes.in.gene.sets.of.interest.rds"))
pathway_annos     <- read_csv(here("data", "input", "my_selected_gene_sets_mapped.csv"),
                              col_select = (c(database, my_name, id)))

TF_annos <- read_tsv(here("data", "input", "Mus_musculus_TF.tsv"))

# Prep data ------------------------------------------------------------

# Preps data for heatmapping by selecting columns, filtering FC,
# ordering columns, splitting FC and CPM data into separate objects,
# and returning FC and CPM data in a list.
prep_hm <- function(data) {

  filt_data <- data |> filter(if_all(contains("Control"), ~ .x > 0.5))
  FC_data <- filt_data |>
    select(contains(c("Name", "Log2FC"))) |>
    relocate(contains("TotalRNA"), .after = "Name") |>
    column_to_rownames("Name")

  CPM_data <- filt_data  |>
    select(contains(c("Name", "Control"))) |>
    relocate(contains("TotalRNA"), .after = "Name") |>
    column_to_rownames("Name")

  output <- list(FCs = FC_data, control_CPMs = CPM_data)
  return(output)
}

RNAdata_deps_hm <- RNAdata_deps |> prep_hm()
RNAdata_deps_DE_hm <- RNAdata_deps_DE |> prep_hm()

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
RNAdata_deps_hm$clust <- RNAdata_deps_hm |> cluster_FC_data(8)
RNAdata_deps_DE_hm$clust <- RNAdata_deps_DE_hm |> cluster_FC_data(8)


# Make global Heatmap -----------------------------------------------------


createDE_Heatmap <- function(data) {

  fc_matrix <- as.matrix(data$FCs)
  cpm_matrix <- as.matrix(data$control_CPMs)
  dep_labels <- colnames(fc_matrix) |> str_split_i("\\.", 1) |> str_replace("^No", "No ")
  cpm_labels <- paste(str_split_i(colnames(cpm_matrix), "\\.", 1),
                      str_split_i(colnames(cpm_matrix), "\\.", 2))
  colnames(cpm_matrix) <- cpm_labels

  k <- data$clust$k

  RNA_type    <- c(rep("Input", 5), rep("Polysome", 5))
  cond_colors <- rep(c("dodgerblue", "red2", "springgreen3", "gold3","purple1"), 2)

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
    col         = list(Expression      = colors_CPM,
                       Cluster = colors_cluster
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
  selected_genes$C8 <- c("Hmgcs1", "Dhcr24", "Sc5d")


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
                     column_split       = RNA_type,
                     row_split = splits,
                     row_title_rot = 45,
                     column_labels = dep_labels,
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
  ht_opt(RESET =  TRUE)
  return(plot)
  }


main_heatmap <- createDE_Heatmap(RNAdata_deps_DE_hm)

saveRDS(main_heatmap, here("output", "objects", "s03_main_heatmap.rds"))

pdf(here("output", "plots", "s03_main_heatmap.pdf"), width = 5, height = 9)
plot(main_heatmap)
    dev.off()


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

    gost_plots <- list()
for (b in names(RNAdata_deps_DE_hm$clust$genes_by_ord_clust)) {
  name <- b
  gost_plots[[name]] <- plot_gostplot(RNAdata_deps_DE_hm$clust$genes_by_ord_clust[b])
}

saveRDS(gost_plots, here("output", "objects", "s03_gost_plots.rds"))
saveRDS(RNAdata_deps_DE_hm$clust$genes_by_ord_clust, here("output", "objects", "s03_clust_genes.rds"))



# Cluster-specific heatmaps -----------------------------------------------

# Update Montana's gene names object with pathway name

pathway_names <- pathway_annos$my_name[which(pathway_annos$id %in% names(pathway_genes))]

names(pathway_genes) <- pathway_names
pathway_genes_tr <- list_transpose(pathway_genes)

createDE_pathway_Heatmap <- function(data, name, folder) {

  cat(name)
  cat("\n")
  cat(nrow(data$FCs))
  cat("\n")

  if (nrow(data$FCs) == 0) {
    message(paste("No detected genes to plot"))
    return(NULL)
  }

  fc_matrix <- as.matrix(data$FCs)
  cpm_matrix <- as.matrix(data$control_CPMs)
  dep_labels <- colnames(fc_matrix) |> str_split_i("\\.", 1) |> str_replace("^No", "No ")
  cpm_labels <- paste(str_split_i(colnames(cpm_matrix), "\\.", 1),
                      str_split_i(colnames(cpm_matrix), "\\.", 2))
  colnames(cpm_matrix) <- cpm_labels


  RNA_type    <- c(rep("Input", 5), rep("Polysome", 5))
  cond_colors <- rep(c("dodgerblue", "red2", "springgreen3", "gold3","purple1"), 2)

  #Main data color mapping
  colors_main    <- colorRamp2(c(-2, 0, 2), c("dodgerblue2", "white", "red1"))
  #Right anno color mapping
  colors_CPM <- colorRamp2(c(0, 4, 8, 15), c("white", "gold1", "orange", "red1"))
  #Build final row annotations
  left_annos <- rowAnnotation(
    Expression         = cpm_matrix,
    #Color assignments
    col         = list(Expression      = colors_CPM
    ),
    border = c(Expression = TRUE),
    annotation_name_offset = unit(2, "mm"),
    annotation_name_gp = gpar(fontsize = 10),
    annotation_width = unit(5, "mm"), width = unit(10, "mm"),
    simple_anno_size_adjust = TRUE,
    show_annotation_name = c(Cluster = FALSE)
  )

  ht_opt$TITLE_PADDING = unit(4, "mm")
  ht_opt$DIMNAME_PADDING = unit(2, "mm")
  ht_opt$ROW_ANNO_PADDING = unit(2, "mm")
  heatmap <- Heatmap(matrix            = fc_matrix,
                     col                = colors_main,
                     left_annotation    = left_annos,
                     name = str_wrap("log2FC vs Control", width = 10),

                     width      = unit(5, "in"),
                     height     = unit((nrow(fc_matrix) / 10), "in"),

                     column_order       = colnames(data$FCs),
                     column_split       = RNA_type,
                     row_title_rot = 45,
                     column_labels = dep_labels,
                     show_row_dend      = FALSE,
                     show_row_names     = TRUE,
                     row_gap            = unit(1, "mm"),
                     column_gap         = unit (2, "mm"),
                     column_names_gp    = gpar(col = cond_colors, fontface = "bold"),
                     row_names_gp       = gpar(fontsize = 8),
                     border             = TRUE,

                     cluster_columns    = FALSE,
                     cluster_row_slices = FALSE,
                     use_raster         = TRUE
  )


  plot <- draw(heatmap, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", column_title = name, column_title_gp = gpar(fontsize=16))

  pdf(here("output", "plots", folder, paste0(name,".pdf")), width = 7, height = (nrow(fc_matrix) / 10) + 5)
  plot(plot)
  dev.off()

  cat("size", (nrow(fc_matrix) / 10))
  cat("\n")

  ht_opt(RESET =  TRUE)
  return(plot)
}




filter_pathway <- function(data, geneset) {
  filt_data <- filter(data, Name %in% geneset)
  return(filt_data)
}


map(pathway_genes_tr$external_gene_name, filter_pathway, data = RNAdata_deps_DE) |>
  map(prep_hm) |>
  walk2(.y = pathway_names, createDE_pathway_Heatmap, folder = "s03_pathway_heatmaps_DE")

all_pdfs <- dir(here("output", "plots", "s03_pathway_heatmaps_DE"), full.names = TRUE)

qpdf::pdf_combine(input = all_pdfs,
                  output = here("output", "plots", "s03_pathway_heatmaps_DE_combined.pdf")
)

map(pathway_genes_tr$external_gene_name, filter_pathway, data = RNAdata_deps) |>
  map(prep_hm) |>
  walk2(.y = pathway_names, createDE_pathway_Heatmap, folder = "s03_pathway_heatmaps_all")

all_pdfs <- dir(here("output", "plots", "s03_pathway_heatmaps_all"), full.names = TRUE)

qpdf::pdf_combine(input = all_pdfs,
                  output = here("output", "plots", "s03_pathway_heatmaps_all_combined.pdf")
)

# Transcription factor heatmaps -------------------------------------------

TF_families <- unique(TF_annos$Family)

TF_genes_by_family <- vector("list", length(TF_families))
#names(TF_genes_by_family) <- TF_families

TF_genes_by_family <- TF_families |>
  set_names()|> # using a named vector here allows the output list to be named too.
  map(\(family) {
    symbol_vec <- TF_annos$Symbol[which(TF_annos$Family == family)]
    return(symbol_vec)
  }
  )



map(TF_genes_by_family, filter_pathway, data = RNAdata_deps_DE) |>
  map(prep_hm) |>
  walk2(.y = TF_families, createDE_pathway_Heatmap, folder = "s03_tf_heatmaps_DE")

TF_pdfs <- dir(here("output", "plots", "s03_tf_heatmaps_DE"), full.names = TRUE)

qpdf::pdf_combine(input = TF_pdfs,
                  output = here("output", "plots", "s03_tf_heatmaps_DE_combined.pdf")
)

map(TF_genes_by_family, filter_pathway, data = RNAdata_deps) |>
  map(prep_hm) |>
  walk2(.y = TF_families, createDE_pathway_Heatmap, folder = "s03_tf_heatmaps_all")

TF_pdfs <- dir(here("output", "plots", "s03_tf_heatmaps_all"), full.names = TRUE)

qpdf::pdf_combine(input = TF_pdfs,
                  output = here("output", "plots", "s03_tf_heatmaps_all_combined.pdf")
)

# Mini pathway-specific and TF heatmaps -----------------------------------

my_tfs <- list()
my_tfs$atf <- c("Atf3", "Atf4", "Atf5", "Atf6")
my_tfs$cebp <- c("Cebpb", "Cebpg", "Ddit3")
my_tfs$ap1 <- c("Fos", "Fosl1", "Fosl2", "Jun", "Junb")
my_tfs$nrf <- c("Nfe2l1", "Nfe2l2")
my_tfs$egr <- c("Egr2", "Egr3")
my_tfs$nr <- c("Nr4a1", "Nr4a3")
my_tfs$irf <- c("Irf1", "Irf7", "Irf9")
my_tfs$stat <- c("Stat1", "Stat3")
my_tfs$nfkb  <- c("Nfkb1", "Nfkb2", "Rel", "Rela")

tfs <- unlist(my_tfs, use.names = FALSE)
if (all(tfs %in% RNAdata_deps$Name) == FALSE) {
  nf <- tfs[!(tfs %in% RNAdata_deps$Name)]
    message(glue("Genes not found: {nf}"))

} else {
  cat("All genes found!")
  }


TF_names <- factor(c("ATF", "C/EBP", "AP-1", "NRF", "EGR", "NR", "IRF", "STAT", "NF-kB"),
                   levels = c("ATF", "C/EBP", "AP-1", "NRF", "EGR", "NR", "IRF", "STAT", "NF-kB"))
TF_splits <- rep(TF_names, lengths(my_tfs))

create_sub_Heatmap <- function(data, name, show_left_anno = TRUE, row_order, splits) {

  fc_matrix <- as.matrix(data$FCs)
  cpm_matrix <- as.matrix(data$control_CPMs)
  dep_labels <- colnames(fc_matrix) |> str_split_i("\\.", 1) |> str_replace("^No", "No ")
  cpm_labels <- paste(str_split_i(colnames(cpm_matrix), "\\.", 1),
                      str_split_i(colnames(cpm_matrix), "\\.", 2))
  colnames(cpm_matrix) <- cpm_labels


  RNA_type    <- c(rep("Input", 5), rep("Polysome", 5))
  cond_colors <- rep(c("dodgerblue", "red2", "springgreen3", "gold3","purple1"), 2)

  #Main data color mapping
  colors_main    <- colorRamp2(c(-2, 0, 2), c("dodgerblue2", "white", "red1"))
  #Right anno color mapping
  colors_CPM <- colorRamp2(c(0, 4, 8, 15), c("white", "gold1", "orange", "red1"))

  if (show_left_anno == TRUE) {  #Build final row annotations
  left_annos <- rowAnnotation(
    Expression         = cpm_matrix,
    #Color assignments
    col         = list(Expression      = colors_CPM
    ),
    border = c(Expression = TRUE),
    annotation_name_offset = unit(2, "mm"),
    annotation_name_gp = gpar(fontsize = 10),
    annotation_width = unit(5, "mm"), width = unit(10, "mm"),
    simple_anno_size_adjust = TRUE,
    show_annotation_name = c(Cluster = FALSE)
  )
  } else {left_annos <- NULL}

  ht_opt$TITLE_PADDING = unit(3, "mm")
  ht_opt$DIMNAME_PADDING = unit(2, "mm")
  ht_opt$ROW_ANNO_PADDING = unit(2, "mm")
  heatmap <- Heatmap(matrix            = fc_matrix,
                     col                = colors_main,
                     left_annotation    = left_annos,
                     name = str_wrap("log2FC vs Control", width = 10),

                     width      = unit(1.8, "in"),
                     height     = unit((nrow(fc_matrix) / 6), "in"),

                     column_order       = colnames(data$FCs),
                     row_order       = row_order,
                     column_split       = RNA_type,
                     row_split       = splits,
                     row_title_rot = 0,
                     column_labels = dep_labels,
                     show_row_dend      = FALSE,
                     show_row_names     = TRUE,
                     row_gap            = unit(0, "mm"),
                     column_gap         = unit (2, "mm"),
                     column_names_gp    = gpar(col = cond_colors,
                                               fontsize = 11.5,
                                               fontface = "bold"),
                     row_names_gp       = gpar(fontsize = 10
                                               ),
                     border             = TRUE,
                     show_heatmap_legend = FALSE,

                     cluster_columns    = FALSE,
                     cluster_rows    = FALSE,
                     cluster_row_slices = FALSE,
                     use_raster         = FALSE
  )


  plot <- draw(heatmap, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", column_title = name, column_title_gp = gpar(fontsize=16))

  pdf(here("output", "plots", paste0(name,".pdf")), width = 4, height = (nrow(fc_matrix) / 6) + 4)
  plot(plot)
  dev.off()
  ht_opt(RESET =  TRUE)
  return(plot)

}

my_TF_heatmap <- RNAdata_deps |>
  filter_pathway(geneset = tfs) |>
  arrange(match(Name, tfs)) |>
  prep_hm() |>
  create_sub_Heatmap(name = "s06_selected_tfs", show_left_anno = FALSE, row_order = tfs, splits = TF_splits)



# Selected pathway heatmaps -----------------------------------------------

genes_up <- list()
genes_up$aa <- c("Gpt2", "Asns", "Got1", "Got2", "Glud1", "Pycr1", "Aldh18a1", "Cth",
       "Mtr", "Bcat1", "Phgdh", "Psat1", "Psph", "Shmt2", "Mthfr", "Mthfd1l", "Mthfd2", "Mthfd1", "Aldh1l2")
genes_up$transport <- c("Slc3a2", "Slc7a1", "Slc7a5", "Slc6a9", "Slc38a1", "Slc38a2", "Slc1a5", "Slc7a11", "Slc1a4")
genes_up$ars <- c("Aars1", "Cars1", "Lars1", "Eprs1", "Nars1", "Vars1", "Mars1")
genes_up$glyc <- c("Hif1a", "Hk1", "Hk2", "G6pdx", "Pcx", "Pck2", "Pkm", "Pfkp", "Tkt", "Shpk", "Pgd", "Eno1", "Eno1b", "Eno2")
genes_up$stress <- c("Gadd45a", "Gadd45b", "Gadd45g", "Chek1", "Fas", "Atm", "Atr")

genes_up_all <- unlist(genes_up, use.names = FALSE)
if (all(genes_up_all %in% RNAdata_deps$Name) == FALSE) {
  nf <- genes_up_all[!(genes_up_all %in% RNAdata_deps$Name)]
  message(glue("Genes not found: {nf}"))

} else {
  cat("All genes found!")
}
length(genes_up_all)
length(unique(genes_up_all))

genes_up_all[duplicated(genes_up_all)]



genes_up_names <- factor(c("AAs", "Transport", "ARSs", "Glc. / PPP", "DNA Stress"),
                   levels = c("AAs", "Transport", "ARSs", "Glc. / PPP", "DNA Stress"))
genes_up_splits <- rep(genes_up_names, lengths(genes_up))


genes_up_heatmap <- RNAdata_deps |>
  filter_pathway(geneset = genes_up_all) |>
  arrange(match(Name, genes_up_all)) |>
  prep_hm() |>
  create_sub_Heatmap(name = "s06_selected_genes_up", show_left_anno = FALSE, row_order = genes_up_all, splits = genes_up_splits)

genes_dn <- list()
genes_dn$cribo <- c( "Rpl36", "Rps16", "Rpl34", "Rpl28", "Rpl35", "Rps7", "Rpl29", "Rps28", "Rps9", "Rps21", "Rpl13", "Rpl27",
                     "Rpl15", "Rpl22", "Rpl31", "Rpl38", "Rpl36a", "Rpl37", "Rps25", "Rps17", "Rpl37a", "Rpl32", "Rpl39")
genes_dn$mribo <- c("Mrpl30", "Mrpl17", "Mrps16", "Mrpl15", "Mrpl14", "Mrpl20", "Mrpl34", "Mrpl36", "Mrps21", "Mrps17", "Mrps14",
                    "Mrpl22", "Mrps10", "Mrps7", "Mrpl16", "Mrps15", "Mrpl27", "Mrpl33", "Mrps18c", "Mrpl28", "Mrpl11", "Mrpl18")
genes_dn$oxphos <- c("Ndufb9", "Ndufa7", "Ndufb1", "Ndufa3", "Ndufa12", "Uqcrfs1", "Uqcrb", "Ndufb3", "Ndufs6",
                     "Ndufab1","Ndufa5","Ndufb10","Cycs","Ndufa4","Ndufb4","Uqcr11","Uqcrh")
genes_dn$chol <- c("Hsd17b7", "Dhcr24", "Sqle", "Fdft1", "Nsdhl", "Sc5d", "Mvk", "Pmvk", "Hmgcs1")

genes_dn_all <- unlist(genes_dn, use.names = FALSE)
if (all(genes_dn_all %in% RNAdata_deps$Name) == FALSE) {
  nf <- genes_dn_all[!(genes_dn_all %in% RNAdata_deps$Name)]
  message(glue("Genes not found: {nf}"))

} else {
  cat("All genes found!")
}
length(genes_dn_all)
length(unique(genes_dn_all))


genes_dn_names <- factor(c("cRibo", "mRibo", "OXPHOS", "Sterols"),
                         levels = c("cRibo", "mRibo", "OXPHOS", "Sterols"))
genes_dn_splits <- rep(genes_dn_names, lengths(genes_dn))

genes_dn_heatmap <- RNAdata_deps |>
  filter_pathway(geneset = genes_dn_all) |>
  arrange(match(Name, genes_dn_all)) |>
  prep_hm() |>
  create_sub_Heatmap(name = "s06_selected_genes_dn", show_left_anno = FALSE, row_order = genes_dn_all, splits = genes_dn_splits)
