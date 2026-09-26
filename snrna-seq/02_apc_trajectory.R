## 02_apc_trajectory.R ---------------------------------------------------------
## Adipocyte progenitor and adipocyte subclustering, and pseudotime ordering.
## This is the analysis behind the main figure; it continues from the object
## saved by 01_clustering.R.
## -----------------------------------------------------------------------------

## ---- packages ---------------------------------------------------------------
## monocle3 and SeuratWrappers are not on CRAN; install them once with:
##   install.packages("remotes")
##   remotes::install_github("cole-trapnell-lab/monocle3")
##   remotes::install_github("satijalab/seurat-wrappers")

for (pkg in c("Seurat", "patchwork", "ggplot2", "dplyr", "openxlsx")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("installing ", pkg, " ...")
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)
library(openxlsx)
library(monocle3)
library(SeuratWrappers)


## ---- settings ---------------------------------------------------------------

IN_RDS  <- "Merge.rds"
OUT_RDS <- "APC_Ad.rds"

## Clusters of the merged object carried forward here: cluster 0 is the Pdgfra+
## progenitor compartment, clusters 1 and 6 the Adipoq+ adipocytes.
PARENT_CLUSTERS <- c("0", "1", "6")

SAMPLE_LEVELS <- c("VcF1_P0", "dVcF1_P0", "VcF1_P7", "dVcF1_P7")

## subclustering
N_DIMS     <- 10
RESOLUTION <- 0.2

## Trajectory roots. Both refer to subclusters of this object and are ordered
## separately, because the two trajectories answer different questions.
ROOTS <- c(committed_progenitors = "0",   # main trajectory
           Dpp4_progenitors      = "2")   # alternative root, upstream population

MARKERS <- c("Dpp4", "Dpp6", "Dpp10", "Esrrg")


## ---- subset and recluster ---------------------------------------------------

merged <- readRDS(IN_RDS)

apc <- subset(merged, cells = WhichCells(
  merged, expression = seurat_clusters %in% PARENT_CLUSTERS))

table(apc$seurat_clusters)

apc <- NormalizeData(apc)
apc <- FindVariableFeatures(apc)
apc <- ScaleData(apc)
apc <- RunPCA(apc, verbose = FALSE)

ElbowPlot(apc, ndims = 30, reduction = "pca")

apc <- FindNeighbors(apc, dims = 1:N_DIMS)
apc <- FindClusters(apc, resolution = RESOLUTION, verbose = FALSE)
apc <- RunUMAP(apc, dims = 1:N_DIMS)
apc <- RunTSNE(apc, dims = 1:N_DIMS)

apc$sample <- factor(apc$sample, levels = SAMPLE_LEVELS)

ggsave("apc_umap.pdf",
       DimPlot(apc, label = TRUE), width = 10, height = 8)
ggsave("apc_umap_by_sample.pdf",
       DimPlot(apc, split.by = "sample", label = TRUE), width = 24, height = 6)
ggsave("apc_tsne_by_sample.pdf",
       DimPlot(apc, reduction = "tsne", split.by = "sample", label = TRUE),
       width = 24, height = 6)


## ---- marker expression ------------------------------------------------------

blank_axes <- theme(axis.title = element_blank(),
                    axis.text  = element_blank(),
                    axis.ticks = element_blank(),
                    axis.line  = element_blank(),
                    panel.border    = element_blank(),
                    legend.position = "none")

ggsave("apc_featureplot.pdf",
       FeaturePlot(apc, features = MARKERS, ncol = length(MARKERS),
                   split.by = "sample") & blank_axes,
       width = 14, height = 12)


## ---- subcluster composition per sample --------------------------------------

Idents(apc) <- "seurat_clusters"

counts  <- as.data.frame.matrix(table(Idents(apc), apc$sample))
counts  <- counts[, SAMPLE_LEVELS[SAMPLE_LEVELS %in% colnames(counts)],
                  drop = FALSE]
percent <- prop.table(as.matrix(counts), margin = 2) * 100
colnames(percent) <- paste0(colnames(percent), "_%")

composition <- cbind(Cluster = rownames(counts), counts, round(percent, 2))
write.xlsx(composition, "apc_composition.xlsx")
print(composition)


## ---- subcluster markers -----------------------------------------------------

apc_markers <- FindAllMarkers(apc, only.pos = TRUE)
write.xlsx(apc_markers, "apc_markers.xlsx")


## ---- pseudotime (monocle3) --------------------------------------------------
## The root is set to the subcluster identified as the least differentiated
## progenitor population; ordering is otherwise unsupervised.

## The graph is learned once; only the root differs between the two orderings.
cds <- as.cell_data_set(apc)
cds <- cluster_cells(cds)
cds <- learn_graph(cds, use_partition = FALSE)

for (root_name in names(ROOTS)) {

  root_cluster <- ROOTS[[root_name]]
  root_cells   <- colnames(apc)[apc$seurat_clusters == root_cluster]
  cds_rooted   <- order_cells(cds, root_cells = root_cells)

  ggsave(sprintf("apc_pseudotime_%s.pdf", root_name),
         plot_cells(cds_rooted,
                    color_cells_by        = "pseudotime",
                    show_trajectory_graph = TRUE,
                    label_cell_groups     = FALSE,
                    label_leaves          = FALSE,
                    label_branch_points   = FALSE,
                    label_roots           = FALSE),
         width = 8, height = 7)

  pt <- pseudotime(cds_rooted)[colnames(apc)]
  apc[[paste0("pseudotime_", root_name)]] <- pt

  write.xlsx(
    data.frame(cell       = colnames(apc),
               sample     = as.character(apc$sample),
               cluster    = as.character(apc$seurat_clusters),
               pseudotime = pt),
    sprintf("apc_pseudotime_%s.xlsx", root_name))

  cat("pseudotime rooted in cluster", root_cluster,
      sprintf("(%s)\n", root_name))
}


## ---- save -------------------------------------------------------------------

saveRDS(apc, file = OUT_RDS)
cat("saved", OUT_RDS, "\n")
